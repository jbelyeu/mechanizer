genome=https://raw.githubusercontent.com/gogetdata/ggd-recipes/master/genomes/Homo_sapiens/GRCh38/GRCh38.genome

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz
#these are the fields in the file downloaded above, which you can check by looking at the sql commands: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.sql
header="bin	chrom	chromStart	chromEnd	name	score	strand	otherChrom	otherStart	otherEnd	otherSize	uid	posBasesHit	testResult	verdict	chits	ccov	alignfile	alignL	indelN	indelS	alignB	matchB	mismatchB	transitionsB	transversionsB	fracMatch	fracMatchIndel	jcK	k2K"
fixed_header="bin	#chrom	start	end	name	score	strand	otherChrom	otherStart	otherEnd	otherSize	uid	posBasesHit	testResult	verdict	chits	ccov	alignfile	alignL	indelN	indelS	alignB	matchB	mismatchB	transitionsB	transversionsB	fracMatch	fracMatchIndel	jcK	k2K"
echo "$fixed_header" | cut -f 2-13,18- > genomicSuperDups.bed

zgrep "." genomicSuperDups.txt.gz \
    | cut -f 2-13,18-\
    | sed 's/chrM/MT/g'\
    | grep -v "chrUn"\
    | grep -v "random"\
    | sed 's/chr//g'\
    | awk '{if (($22 > 0.95) && ($1 == $7) && (($8-$3) > 100)) {print $0}}'\
    | gsort /dev/stdin $genome\
    >>genomicSuperDups.bed
bgzip genomicSuperDups.bed
tabix genomicSuperDups.bed.gz

rm genomicSuperDups.txt.gz
