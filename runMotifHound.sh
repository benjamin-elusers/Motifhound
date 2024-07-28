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
	X=$((((L-2)/2)+1));
	echo "============"
	echo "L=$L X=$X n=$n N=$N "
	echo "============"
	
	set=$(basename $S)
	back=$(basename $B)
	motif=${OUT}/${set##.fasta}.L${L}_X$X.motif 
	scanned=${OUT}/${back##.fasta}.L${L}_X$X.scanned 
	enriched=${OUT}/${set##.motif}.L${L}_X$X.enriched

	if [ -f $motif ]; then
		echo "SKIPPED: ALREADY ENUMERATED MOTIFS"
	else
		echo "$BIN/MotifEnumeration.exe -f $S -o $motif -m $L -k 5  -a 'NT' -x 0 -y $X "
		$BIN/MotifEnumeration.exe -f $S -o $motif -m $L -k 5  -a 'NT' -x 0 -y $X	
	fi


	if [ -f $scanned ];
		echo "SKIPPED: ALREADY SCANNED MOTIFS"
	elif [ ! -f $scanned ] && [ -f $motif ]; then
		echo "$BIN/MotifScan.exe -e $motif -f $B -o $scanned -a 'NT' "
		$BIN/MotifScan.exe -e $motif -f $B -o $scanned -a 'NT'
	else 
		echo "ERROR: REQUIRE MOTIF FILE IN SET"
	fi
	
	if [ -f $enriched ]; then
		echo "SKIPPED: ALREADY COMPUTED MOTIF ENRICHMENT"
	elif [ ! -f $enriched ] && [ -f $scanned ] && [ -f $motif ]; then
		echo "$BIN/MotifEnrichment.exe -s $motif -b $scanned -o $enriched -n $n -N $N"
		$BIN/MotifEnrichment.exe -s $motif -b $scanned -o $enriched -n $n -N $N
	else
		echo "ERROR: REQUIRE MOTIF FILES IN BOTH SET AND BACKGROUND"
	fi
done