#! /usr/bin/env python
from __future__ import print_function
# Python 2/3 compatibility
import sys
import argparse
from cyvcf2 import VCF
import pandas as pd
import pysam
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
MATCH = 1
MISMATCH = -0.5
GAP_OPEN = -0.25
GAP_EXTEND = -0.05

def setup_args():
    parser = argparse.ArgumentParser(description="guesses the most likely mutational mechanism for structural variants in a VCF")
    parser.add_argument(
        "-b",
        "--bed",
        help="bed of structural variants with header that labels the chrom,start,end, and svtype columns",
        required=False,
        default="/Users/jon/Research/scripts/de_novo_sv/ceph_denovos/test/data/dnms.bed"
    )
    parser.add_argument(
        "-v", 
        "--vcf", 
        help="vcf of structural variants",
        required=False
    )
    parser.add_argument(
        "-r", 
        "--ref",
        help="reference genome",
    #    required=True,
        default="/Users/jon/miniconda3/share/ggd/Homo_sapiens/GRCh37/grch37-reference-genome-ensembl-v1/1/grch37-reference-genome-ensembl-v1.fa"
    )
    args = parser.parse_args()

    if not args.vcf and not args.bed:
        sys.exit("ERROR: mechanizer requires either --bed or --vcf argument")
    elif None not in [args.vcf, args.bed]:
        sys.exit("ERROR: mechanizer requires either --bed or --vcf argument, not both")
    return args


def main(args):
    svs = []
    header = []
    with open(args.bed, 'r') as sv_file:
        for line in sv_file:
            fields = line.strip().split()
            if line[0] == "#":
                if "svtype" in line.lower():
                    header = [f.lower().strip("#") for f in fields]
                continue
            svs.append(fields)

    try:
        chrom_idx = header.index("chrom")
        start_idx = header.index("start")
        end_idx = header.index("end")
        svtype_idx = header.index("svtype")
    except Exception as e:
        sys.exit("header missing required field. Must contain chrom, start, end, svtype")

    reference = pysam.FastaFile(args.ref)
    homology_sizes = [3, 10, 100, 1000]
    print("\t".join([str(s) for s in homology_sizes]))
    for sv in svs:
        chrom = sv[chrom_idx]
        start = int(sv[start_idx])
        end = int(sv[end_idx])
        svtype = sv[svtype_idx]
        matches = []
        #check for homologies of different sizes
        for match_size in homology_sizes:
            upstream = reference.fetch(chrom, start-match_size, start)
            downstream = reference.fetch(chrom, end, end+match_size)

            align_score = pairwise2.align.globalms(upstream, downstream, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND, score_only=True)
            matches.append(align_score/match_size)
            
        print("\t".join([str(m) for m in matches]))


if __name__ == "__main__":
    args = setup_args()
    main(args)


