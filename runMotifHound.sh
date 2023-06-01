#! /bin/sh

CURRENT=$(dirname $0);
BIN="./bin"
DAT="data"
OUT="output"
SRC="src/V2/"

gcc $SRC/MotifEnumeration.c $SRC/safealloc.c $SRC/fasta.c $SRC/motif.c -o $BIN/MotifEnumeration.exe -lm -lJudy -lpthread
gcc $SRC/MotifScan.c        $SRC/safealloc.c $SRC/fasta.c $SRC/motif.c -o $BIN/MotifScan.exe        -lm -lJudy -lpthread
gcc $SRC/MotifEnrichment.c  $SRC/safealloc.c $SRC/fasta.c $SRC/motif.c -o $BIN/MotifEnrichment.exe  -lm -lJudy -lpthread


S="${DAT}/COCO_pairs.fasta"
n=$(grep -c ">" $S)
B="${DAT}/POST_pairs.fasta"
N=$(grep -c ">" $B)

for L in 8 12 16 20; do
	X=$(((L-2)/2));
	echo $L $X $n $N;
	set=$(basename $S)
	back=$(basename $B)
	motif=${OUT}/${set##.fasta}.L${L}_X$X.motif 
	scanned=${OUT}/${back##.fasta}.L${L}_X$X.scanned 
	enriched=${OUT}/${set##.motif}.L${L}_X$X.enriched
	echo "$BIN/MotifEnumeration.exe -f $S -o $motif -m $L -k 3  -a 'NT' -x 0 -y $X "
	#$BIN/MotifEnumeration.exe -f $S -o $motif -m $L -k 3  -a 'NT' -x 0 -y $X
	
	echo "$BIN/MotifScan.exe -e $motif -f $B -o $scanned -a 'NT' "
	#$BIN/MotifScan.exe -e $motif -f $B -o $scanned -a 'NT'
	
	echo "$BIN/MotifEnrichment.exe -s $motif -b $scanned -o $enriched -n $n -N $N"
	$BIN/MotifEnrichment.exe -s $motif -b $scanned -o $enriched -n $n -N $N
done