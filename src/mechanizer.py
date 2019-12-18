#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys
import argparse
import pandas as pd
import pysam
import gzip
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
MATCH = 1
MISMATCH = -0.5
GAP_OPEN = -0.25
GAP_EXTEND = -0.05
SPLITTER_ERR_MARGIN=3
MIN_MAPQ=10
HOMOLOGY_CHECK_DIST = 100

def setup_args():
    parser = argparse.ArgumentParser(description="guesses the most likely mutational mechanism for structural variants in a VCF")
    parser.add_argument(
        "-b",
        "--bed",
        help="bed of structural variants with header that labels the chrom,start,end, and svtype columns",
        required=True,
    )
    parser.add_argument(
        "-r", 
        "--ref",
        help="reference genome",
        required=True,
    )
    parser.add_argument(
        "-p", 
        "--project_name",
        help="name of the proejct",
        required=True,
    )
    parser.add_argument(
        "--repeats",
        help="bed file of repeats that might act as NAHR substrates, where the header of the bed file identifies the coordinates of `chrom,start,end,otherChrom,otherStart,otherEnd`",
        required=True,
    )
    args = parser.parse_args()

    if not args.bed:
        sys.exit("ERROR: mechanizer requires either --bed argument")
    return args

def measure_homology(ref_file, chrom, start, end):
    reference = pysam.FastaFile(ref_file)
    
    #matches = []
    #check for homologies of different sizes
    #for match_size in homology_sizes:
    #instead of looking for a fuzzy match, let's find a hard one by counting the bases to the first mismatch
    upstream = None
    downstream = None
    chrom = str(chrom)
    try:
        upstream = reference.fetch(chrom, start-HOMOLOGY_CHECK_DIST, start)[HOMOLOGY_CHECK_DIST::-1]
        downstream = reference.fetch(chrom, end, end+HOMOLOGY_CHECK_DIST)
    except KeyError:
        if "chr" in chrom:
            chrom = chrom.strip("chr")
        else:
            chrom = "chr"+chrom
        upstream = reference.fetch(chrom, start-HOMOLOGY_CHECK_DIST, start)[HOMOLOGY_CHECK_DIST::-1]
        downstream = reference.fetch(chrom, end, end+HOMOLOGY_CHECK_DIST)
    matchlen = 0

    for i,downbase in enumerate(downstream):
        if downbase == upstream[i]:
            matchlen += 1
        else:
            break
    return matchlen


def goodread(read):
    if not read:
        return False
    if (read.is_qcfail
            or read.is_unmapped
            or read.is_duplicate
            or int(read.mapping_quality) < MIN_MAPQ
            or read.is_secondary
            or read.is_supplementary
            or (read.next_reference_id != read.reference_id)
    ):
        return False
    return True

def get_reads(bam_filehandle, chrom, start, end):
    bamit = None
    try:
        bamit = bam_filehandle.fetch(chrom, start, end)
    except KeyError:
        if "chr" in chrom:
            chrom = chrom.strip("chr")
        else:
            chrom = "chr"+chrom
        bamit = bam_open.fetch(chrom, start, end)
    return [read for read in bamit]



def get_splits(bamfile, chrom, start, end):
    splits = []
    bam_filehandle = pysam.AlignmentFile(bamfile)
    chrom = str(chrom)
    reads = get_reads(bam_filehandle, chrom, max(0,start-SPLITTER_ERR_MARGIN), start+SPLITTER_ERR_MARGIN)
    reads += get_reads(bam_filehandle, chrom, max(0,end-SPLITTER_ERR_MARGIN), end+SPLITTER_ERR_MARGIN)

    for read in reads:
        ref_positions = read.get_reference_positions()
        if len(ref_positions) == 0:
            continue
        right_clip = max(ref_positions)
        left_clip = min(ref_positions)
        if (goodread(read) 
                and read.has_tag("SA") 
                and ((right_clip-SPLITTER_ERR_MARGIN) <= start <= (right_clip+SPLITTER_ERR_MARGIN)
                    or (left_clip-SPLITTER_ERR_MARGIN) <= start <= (left_clip+SPLITTER_ERR_MARGIN)
                    or (right_clip-SPLITTER_ERR_MARGIN) <= end <= (right_clip+SPLITTER_ERR_MARGIN)
                    or (left_clip-SPLITTER_ERR_MARGIN) <= end <= (left_clip+SPLITTER_ERR_MARGIN))):
            splits.append(read)
    return splits

def get_clipping(read):
    """
    find the largest clipped region in the read and return the coordinates thereof
    """
    positions = read.get_reference_positions(full_length=True)
    #find the longest run of None values in the positions
    curr_start = None
    longest_len = 0
    longest_start = 0
    for i,pos in enumerate(positions):
        #if curr_start isn't set, we haven't started a gap
        if curr_start is None:
            if pos is None: 
                #we need to start a gap
                curr_start = i
            else: 
                #we don't need to do anything
                continue
        else: #we're in a gap
            if pos is None: 
                if curr_start == longest_start: #check if this is already the longest gap
                    #just extend the gap
                    longest_len+=1
                elif (i - curr_start) > longest_len: #not yet set as longest gap but longer than the previous longest
                    longest_len = i - curr_start
                    longest_start = curr_start
                else: #not the longest gap yet, so move on
                    continue
            else: #we're in a gap, but it has ended
                curr_start = None
    
    #get the reference positions of the gap start and end
    clipstart = -1
    if longest_start == 0:
        clipstart = read.reference_start
    
    clipend = -1
    if longest_len+longest_start == 150:
        clipend = read.reference_end

    return clipstart,clipend 

def is_microhomology_mediated(sv, ref):
    left_supporting_splits = 0
    right_supporting_splits = 0
    non_supporting_splits = 0
    splits = get_splits(sv.loc["bam"], sv.loc["chrom"], sv.loc["start"], sv.loc["end"])
    for split in splits:
        clipstart,clipend = get_clipping(split)
        if (((sv.loc["start"] - SPLITTER_ERR_MARGIN) <= clipstart <= (sv.loc["start"]+SPLITTER_ERR_MARGIN)) or
                ((sv.loc["start"] - SPLITTER_ERR_MARGIN) <= clipend <= (sv.loc["start"]+SPLITTER_ERR_MARGIN))):
            left_supporting_splits += 1
        elif (((sv.loc["end"] - SPLITTER_ERR_MARGIN) <= clipstart <= (sv.loc["end"]+SPLITTER_ERR_MARGIN)) or
                ((sv.loc["end"] - SPLITTER_ERR_MARGIN) <= clipend <= (sv.loc["end"]+SPLITTER_ERR_MARGIN))):
            right_supporting_splits += 1
        else:
            non_supporting_splits += 1
    
    microhomology = 0
    if (left_supporting_splits > 2
            and right_supporting_splits > 2 
            and (left_supporting_splits+right_supporting_splits) > non_supporting_splits):
        #we can be quite confident of the breakpoints
        microhomology = measure_homology(args.ref, sv.loc["chrom"], sv.loc["start"], sv.loc["end"])
    return (microhomology>2)



def parse_svs(bed):
    svs = []
    header_idxs = {}
    required_fields = ["chrom","start","end","svtype","kid_id","bam"]
     
    with open(bed, 'r') as sv_file:
        for i,line in enumerate(sv_file):
            if line[0] == "#" and i == 0:
                header_fields = line.lower().strip("#").split()
                for required_field in required_fields:
                    if not required_field in header_fields:
                        sys.exit("missing required header field "+required_field+"in repeats file")
                    header_idxs[required_field] = header_fields.index(required_field)
            elif line[0] != "#":
                sv = []
                fields = line.strip().split()
                for i,required_field in enumerate(required_fields):
                    sv.append(to_int(fields[header_idxs[required_field]]))
                svs.append(sv)
    return pd.DataFrame(svs,columns=required_fields)



def to_int(field):
    try:
        return int(field)
    except ValueError:
        return field



def parse_repeats(bed):
    encoding = "utf-8"
    repeats = []
    required_fields = ['chrom','start','end','otherchrom','otherstart','otherend']
    with gzip.open(bed, 'r') as repeat_file:
        header_idxs = {}
        for line in repeat_file:
            line = line.decode(encoding)
            if line[0] == "#":
                header_fields = line.lower().strip("#").split()
                for required_field in required_fields:
                    if not required_field in header_fields:
                        sys.exit("missing required header field "+required_field+"in repeats file")
                    header_idxs[required_field] = header_fields.index(required_field)
            else:
                repeat_entry = []
                fields = line.strip().split()
                for i,required_field in enumerate(required_fields):
                    repeat_entry.append(to_int(fields[header_idxs[required_field]]))
                repeats.append(repeat_entry)
    return pd.DataFrame(repeats, columns=required_fields)

def has_flanking_repeats(sv, repeat_file):
    if sv.loc['svtype'] not in ["DEL","DUP"]:
        return False

    repeats = parse_repeats(repeat_file)
    # find all the repeat pairs where 
    # the first of the pair ends within 100 bp before the SV start and not after SV start
    # and the second of the pair starts after the SV end and within 100 bp after the SV end
    search_dist = int((sv.loc['end']-sv.loc['start'])/20)

    repeats["chrom"] = repeats['chrom'].astype(str)
    candidates = repeats[
              ((repeats["chrom"] == str(sv.loc['chrom']).strip("chr")) 
                  | (repeats["chrom"] == "chr"+str(sv.loc['chrom']))) 
            & (repeats['end'] >= (sv.loc['start']-search_dist))
            & (repeats['end'] <= sv.loc['start']+search_dist)
            & (repeats['otherstart'] >= sv.loc['end']-search_dist)
            & (repeats['otherstart'] <= (sv.loc['end']+search_dist))
    ]

    return (len(candidates) > 0)



def main(args):
    svs = parse_svs(args.bed)
    svs = svs[(svs['end']-svs['start']) > 50]
    likely_mechanisms = []
        
    for idx,sv in svs.iterrows():
        if (sv.loc['end']-sv.loc['start']) < 50:
            continue
        #check for breakpoint microhomology: replication-based errors
        repli_based_likely = is_microhomology_mediated(sv, args.ref)
        
        #check for breakpoint macrohomology: NAHR
        nahr_likely = has_flanking_repeats(sv, args.repeats)
        
        mechanism = ""
        if nahr_likely:
            mechanism = "NAHR"
        elif repli_based_likely:
            mechanism = "RBM"
        else:
            mechanism = "NHEJ"
        likely_mechanisms.append(mechanism)
    svs["pred_mechanism"] = likely_mechanisms
    svs.sort_values(by=['chrom','start','end'],inplace=True)

    svs.rename(columns={'chrom':'#chrom'}, inplace=True)
    with open("{}.bed".format(args.project_name), 'w') as res_file:    
        print(svs.to_csv(header=True, sep="\t", index=False), file=res_file)


if __name__ == "__main__":
    args = setup_args()
    main(args)
