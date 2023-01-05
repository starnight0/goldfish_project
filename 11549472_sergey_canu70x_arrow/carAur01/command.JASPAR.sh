# WD: ../JASPAR_run1
mkdir fimo_split
>r01.fimo.sh
IFS=$'\n'
for line in `cat /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.withM.sm.fa.by_ctg/files`
do
ctg=`echo $line | cut -f 1`
fa=`echo $line | cut -f 2`
echo 'fimo --o fimo_split/'$ctg' ~/data/JASPAR/JASPAR2018/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt "'$fa'"' >> r01.fimo.sh
done
#cmd=r01.fimo; if ! [ -d $cmd.log ]; then mkdir $cmd.log; fi; swarm --time=2:00:00 -g 8 -b 100 --logdir $cmd.log -f $cmd.sh -m meme

mkdir mast_split
>r01.mast.sh
IFS=$'\n'
for line in `cat /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.withM.sm.fa.by_ctg/files`
do
ctg=`echo $line | cut -f 1`
fa=`echo $line | cut -f 2`
echo 'mast --o mast_split/'$ctg' ~/data/JASPAR/JASPAR2018/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt "'$fa'"' >> r01.mast.sh
done
#cmd=r01.mast; if ! [ -d $cmd.log ]; then mkdir $cmd.log; fi; swarm --time=2:00:00 -g 8 -b 100 --logdir $cmd.log -f $cmd.sh -m meme

>r01.mast.x1.sh
IFS=$'\n';
for line in `cat /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.withM.sm.fa.by_ctg/files`
do
ctg=`echo $line | cut -f 1`
fa=`echo $line | cut -f 2`
if ! [ -d mast_split/$ctg ] || ! [ -f mast_split/$ctg/mast.txt ];
then
echo 'mast --o mast_split/'$ctg' ~/data/JASPAR/JASPAR2018/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt "'$fa'"' >> r01.mast.x1.sh
fi
done

mkdir fimo_out
>fimo_out/fimo.txt;
>fimo_out/fimo.gff;
IFS=$'\n';
for line in `cat /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.withM.sm.fa.by_ctg/files`
do
ctg=`echo $line | cut -f 1`
fa=`echo $line | cut -f 2`
if [ -d fimo_split/$ctg ];
then
    cat fimo_split/$ctg/fimo.txt >> fimo_out/fimo.txt
    cat fimo_split/$ctg/fimo.gff >> fimo_out/fimo.gff
else
    echo $ctg Fail.
fi
done
gzip fimo_out/fimo.txt
gzip fimo_out/fimo.gff
#zcat fimo_out/fimo.txt.gz | awk -F$'\t' -v OFS=$'\t' '{chr=$3; strand=$6; b=$4-1; e=$5; score=$7; p=$8; id1=$1; id2=$2; q=$9; seq=$10; if (p<0.001) { if (p<0.000001) {col="0,236,255";} else if (p<0.00001) { col="49,233,255"; } else if (p<0.0001) {col="110,220,255"} else {col="115,204,216";} print chr,b,e,id2,score,strand } }' | LC_COLLATE=C sort -k1,1 -k2,2n > fimo_out/fimo.bed 
zcat fimo_out/fimo.txt.gz | awk -F$'\t' -v OFS=$'\t' '/^[^#]/ {chr=$3; strand=$6; b=$4-1; e=$5; score=$7; p=$8; id1=$1; id2=$2; q=$9; seq=$10; if (p<0.00001) { if (p<0.000001) {col="49,233,255";} else { col="0,236,236"; } print chr,b,e,id2,score,strand,b,e,col } }' | LC_COLLATE=C sort -S 4G -k1,1 -k2,2n > fimo_out/fimo.p5.bed 
zcat fimo_out/fimo.txt.gz | awk -F$'\t' -v OFS=$'\t' '/^[^#]/ {chr=$3; strand=$6; b=$4-1; e=$5; score=$7; p=$8; id1=$1; id2=$2; q=$9; seq=$10; if (p<0.0001) { if (p<0.000001) {col="110,220,255";} else if (p<0.00001) { col="49,233,255"; } else {col="0,236,255"} print chr,b,e,id2,score,strand,b,e,col } }' | LC_COLLATE=C sort -S 4G -k1,1 -k2,2n > fimo_out/fimo.p4.bed 
cat fimo_out/fimo.p5.bed | awk -F$'\t' -v OFS=$'\t' '{ if ($5>1000) $5=1000; else if ($5<0) $5=0; else {$5=int($5+0.5);} print $0}' > fimo_out/fimo.1.bed
liftOver  fimo_out/fimo.1.bed ../carAur03/carAur01_to_carAur03.liftOver.chain  stdout  fimo_out/carAur03.fimo.1.bed.unmap | LC_COLLATE=C sort -k1,1 -k2,2n > fimo_out/carAur03.fimo.1.bed

cat fimo_out/fimo.p4.bed | awk -F$'\t' -v OFS=$'\t' '{ if ($5>1000) $5=1000; else if ($5<0) $5=0; else {$5=int($5+0.5);} print $0}' | liftOver stdin ../carAur03/carAur01_to_carAur03.liftOver.chain  stdout  fimo_out/carAur03.fimo.p4.bed.unmap  | LC_COLLATE=C sort -S 4G -k1,1 -k2,2n > fimo_out/carAur03.fimo.p4.bed

bedToBigBed -type=bed9 -tab fimo_out/fimo.1.bed ../carAur01/carAur01.withM.sm.fa.fai fimo_out/fimo.bb
bedToBigBed -type=bed9 -tab fimo_out/carAur03.fimo.1.bed ../carAur03/carAur03.sm.fa.fai fimo_out/carAur03.fimo.bb
bedToBigBed -type=bed9 -tab fimo_out/carAur03.fimo.p4.bed ../carAur03/carAur03.sm.fa.fai fimo_out/carAur03.fimo.p4.bb
cp fimo_out/carAur03.fimo.bb  ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.JASPAR.fimo.bb 
cp fimo_out/carAur03.fimo.p4.bb  ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.JASPAR.fimo.p4.bb 
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.JASPAR.fimo.*bb 
