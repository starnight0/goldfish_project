cwd=`pwd`
datadir=`echo $cwd | sed 's/data_cmd/data/'`
cd $datadir
echo '##gff-version 3' > goldfish.arrow.renamed.masked.tRNA.gff
genome1fa=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/goldfish.arrow.renamed.masked.fasta
genome2fa=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur02/carAur02.fa
tail -n +4 goldfish.arrow.renamed.masked.tRNA.out \
| awk -v OFS=$'\t' '{ a="ID="$1"_"$2";Name="$5"_"$6";type="$5";anti="$6";intron_begin="$7";intron_end="$8";"; if ($4>$3) {print $1,"tRNAscanSE","tRNA",$3,$4,$9,"+",".",a} else {print $1,"tRNAscanSE","tRNA",$4,$3,$9,"+",".",a} }' >> goldfish.arrow.renamed.masked.tRNA.gff
ln -sf goldfish.arrow.renamed.masked.tRNA.gff carAur01.tRNA.gff
fn=carAur01.tRNA
gff3ToGenePred -geneNameAttr=Name $fn.gff stdout | genePredToBigGenePred stdin stdout | sort -k1,1 -k2,2n > $fn.biggp
bedToBigBed -extraIndex=name2 -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as $fn.biggp $genome1fa.fai $fn.bb 
cd $cwd
