cwd=`pwd`
datadir=`echo $cwd | sed 's/data_cmd/data/'`

module load ucsc bedtools

cd $datadir

asm1=carAur01
asm1_dir=../carAur01
asm2=carAur03
asm2_dir=.
genome2=$asm2
genome2fa=$asm2.fa
kent=/data/genome/jksrc_v352/kent/src/hg/lib/
chain=carAur01_to_carAur03.liftOver.chain
mkdir -p bbi
mkdir -p big

################################
# transform chain to bed
################################
cat $chain | grep '^chain' | awk -v OFS=$'\t' '$10~/^-/ {print $8,$9-$12,$9-$11,$3":"$4":"$6+1"-"$7,"100","-"} $10!~/^-/ {print $8,$11,$12,$3":"$4":"$6+1"-"$7,"100","+"}' | sort -k1,1 -k2,2n > big/carAur03.assembly.bed

###########################
# repeatmasker
# liftOver repeatmasker
###########################
liftOver -multiple $asm1_dir/big/$asm1.repeatmasker.bed $chain stdout $asm2_dir/$asm2.repeatmasker.bed.unmap | sort -k1,1 -k2,2n > $asm2_dir/big/$asm2.repeatmasker.bed
bedtools merge -i big/carAur03.repeatmasker.bed | cut -f 1-3 > big/carAur03.repeatmasker.merged.bed3

###########################
# mask carAur03.fa
###########################
bedtools maskfasta -soft -fi carAur03.fa -bed big/carAur03.repeatmasker.bed -fo carAur03.sm.fa
faSize -detailed carAur03.sm.fa | sort -k1,1 -k2,2n > carAur03.sm.fa.chromSizes
faToTwoBit carAur03.sm.fa carAur03.2bit
samtools faidx carAur03.sm.fa
mkdir carAur03.sm.fa.by_chr
~/my_program3/src/utility/czl_fasta_split -m 5 -i carAur03.sm.fa -o carAur03.sm.fa.by_chr/ -d 1 --cut-name -1
mkdir carAur03.sm.fa.by_chr/chrUn
mv carAur03.sm.fa.by_chr/tig* carAur03.sm.fa.by_chr/chrUn/
cat carAur03.sm.fa.by_chr/chrUn/*.fa > carAur03.sm.fa.by_chr/chrUn.fa
~/my_program3/src/utility/czl_fasta_split -m 2 -n 50000000 -i carAur03.sm.fa.by_chr/chrUn.fa -o carAur03.sm.fa.by_chr/chrUn. --cut-name -1
#----------------------------------
# GC track
# hgGcPercent -wigOut -doGaps -file=stdout -win=200 -verbose=0 carAur02 carAur02.2bit | wigEncode stdin gc.wig gc.wib
hgGcPercent -wigOut -doGaps -file=stdout -win=100 -verbose=0 carAur03 carAur03.2bit | wigToBigWig stdin carAur03.sm.fa.chromSizes bbi/carAur03.gc.bw
#----------------------------------

# copy agp from broken_contig2/carAur03/carAur03.withM.agp
# cp ../broken_contig2/carAur03/carAur03.withM.agp carAur03.agp
# This agp is incorrect because bug in 'Chromonomer'
# chain to agp
chainSwap carAur01_to_carAur03.liftOver.chain stdout | grep '^chain' | sort -t' ' -k3,3 -k6,6n | awk -v OFS=$'\t' -v i=0 '$3~/^LG/ { if ($3==tchr) {i++; } else {i=1; tchr0=""; te0=0;} tchr=$3; tb=$6+1; te=$7; if($10!="-") {qb=$11+1; qe=$12;} else {qb=$9-$12+1; qe=$9-$11;} qchr=$8; if (te0>0) {print tchr0,te0+1,tb-1,i,"U",tb-te0-1,"scaffold","yes","map"; i++;} print tchr,tb,te,i,"W",qchr,qb,qe,$10; tchr0=tchr; te0=te; }' > carAur03.LG.agp


###########################
# liftOver gene track
###########################
CrossMap.py gff $chain ../carAur01/big/carAur01.gene.gff  big/carAur03.noM.gene.gff 
CrossMap.py gff $chain ../carAur01/big/carAur01.noM.gene.unmasked.gff  big/carAur03.noM.gene.unmasked.gff 
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.gene.bgp $chain big/carAur03.noM.gene.bgp big/carAur03.noM.gene.bgp.unmap 
cat big/carAur03.noM.gene.bgp  | sort -k1,1 -k2,2n > a; mv a big/carAur03.noM.gene.bgp 
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.gene.unmasked.bgp $chain stdout big/carAur03.noM.gene.unmasked.bgp.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.gene.unmasked.bgp
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.gene.masked.bgp $chain stdout big/carAur03.noM.gene.masked.bgp.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.gene.masked.bgp
liftOver -tab $asm1_dir/big/$asm1.gene.bed $chain big/carAur03.noM.gene.bed big/carAur03.noM.gene.bed.unmap
liftOver -tab $asm1_dir/big/$asm1.transcript.bed $chain big/carAur03.noM.transcript.bed big/carAur03.noM.transcript.bed.unmap
liftOver -tab $asm1_dir/big/$asm1.exon.bed $chain big/carAur03.noM.exon.bed big/carAur03.noM.exon.bed.unmap
liftOver -tab $asm1_dir/big/$asm1.cds.bed $chain big/carAur03.noM.cds.bed big/carAur03.noM.cds.bed.unmap

liftOver -tab ../carAur01/big/carAur01.gene.gff.bed $chain big/carAur03.noM.gff.bed big/carAur03.noM.gene.gff.bed.unmap

cat big/carAur03.noM.gene.masked.bgp | sort -k18,18 -k1,1 -k2,2n | awk -F$'\t' -v OFS=$'\t' -v g="" '{ if (g==$18) {if (b>$2) {b=$2;} if(e<$3){e=$3;} } else {if (g) {print chr,b,e,g,s,ss}}  chr=$1;b=$2;e=$3;if ($19=="" || $19==".") g=$18; else g=$18"__"$19;s=$5;ss=$6;} END{print chr,b,e,g,s,ss;}' | sort -k1,1 -k2,2n > big/carAur03.noM.gene.masked.bed
cat big/carAur03.noM.gene.unmasked.bgp | sort -k18,18 -k1,1 -k2,2n | awk -F$'\t' -v OFS=$'\t' -v g="" '{ if (g==$18) {if (b>$2) {b=$2;} if(e<$3){e=$3;} } else { if(g) {print chr,b,e,g,s,ss} } chr=$1;b=$2;e=$3;if ($19=="" || $19==".") g=$18; else g=$18"__"$19;s=$5;ss=$6;} END{print chr,b,e,g,s,ss;}' | sort -k1,1 -k2,2n > big/carAur03.noM.gene.unmasked.bed

###########################
# liftOver ncRNA track
###########################
#cat ../carAur01/big/carAur01.ncrna.gp | awk -F$'\t' '$2!~/tig00037772/' | liftOver -genePred stdin $chain big/carAur03.noM.ncrna.gp big/carAur03.noM.ncrna.gp.unmap
cat ../carAur01/big/carAur01.ncrna.bgp | awk -F$'\t' '$1!~/tig00037772/' | liftOver -tab -bedPlus=12 stdin $chain stdout big/carAur03.noM.ncrna.bgp.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.ncrna.bgp 
#cat big/carAur03.noM.ncrna.bgp | perl -ne 'chomp; my @t=split /\t/; my $n=$t[9]; my @sz=split /,/,$t[10]; my $l=0; my $m=0; foreach my $sz (@sz) {$l+=$sz; if ($l<30) {$m++; }}; $l/=$n; if ($l<50 || $m>=2) { next; } else {print $_, "\n"; }'
cat big/carAur03.noM.ncrna.bgp | awk -F$'\t' '$18!~/^CARNA[0-9]+/' >  big/carAur03.noM.ncrna.f.bgp 
cat big/carAur03.noM.ncrna.f.bgp | awk -F$'\t' '{if (chr!=$1 || e<=$2) {print; chr=$1; e=$3;}}' > big/carAur03.noM.ncrna.f.not_ovl.bgp
cat big/carAur03.noM.ncrna.f.not_ovl.bgp | awk 'tolower($18)~/(mir|let|lin)/' > big/carAur03.noM.ncrna.f.not_ovl.mir.bgp 
cat big/carAur03.noM.ncrna.f.not_ovl.bgp | perl -ne 'my @t=split /\t/; chomp;
my ($chr,$b,$e,$tid) = @t[0..3];  my $s = $t[5];
my $n=$t[9];
my @sz = split /,/, $t[10];
my @bs = split /,/, $t[11];
for (my $i=0; $i<$n; $i++) 
{ print join("\t",($chr,$b+$bs[$i],$b+$bs[$i]+$sz[$i],"$tid:".($i+1),".",$s) ), "\n"; }
' | sort -k1,1 -k2,2n > big/carAur03.noM.ncrna.f.not_ovl.exon.bed

###########################
# liftOver tRNA track
###########################
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.tRNA.bgp $chain stdout big/carAur03.tRNA.bgp.unmap | sort -k1,1 -k2,2n > big/carAur03.tRNA.bgp

###########################
# liftOver miRNA track
###########################
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.miRNA.bigpsl $chain stdout big/carAur03.noM.miRNA.bigpsl.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.miRNA.bigpsl
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.miRNA.I95.bigpsl $chain stdout big/carAur03.noM.miRNA.I95.bigpsl.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.miRNA.I95.bigpsl
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.miRNA.f.bed12 $chain stdout big/carAur03.noM.miRNA.bed.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.miRNA.bed12
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.miRNA.I95.f.bed12 $chain stdout big/carAur03.noM.miRNA.I95.bed12.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.miRNA.I95.bed12

###########################
# liftOver Rfam track
###########################
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.rfam.bgp $chain stdout big/carAur03.noM.rfam.bgp.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.rfam.bgp

###########################
# liftOver Interpro track
###########################
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.ips.bgp $chain stdout big/carAur03.noM.ips.bgp.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.ips.bgp
liftOver -tab -bedPlus=12 ../carAur01/big/carAur01.noM.GO.rename.bgp $chain stdout big/carAur03.noM.GO.bgp.unmap | sort -k1,1 -k2,2n > big/carAur03.noM.GO.bgp

###########################
# liftOver est2genome track
###########################
mkdir big/est2genome
mkdir bbi/est2genome
chain=carAur01_to_carAur03.liftOver.chain
for sm in `cat ../Trinity_GG_run01/sample`
do
f=big/est2genome/$sm.trinity_GG
liftOver -bedPlus=12 -tab ../carAur01/$f.bigpsl $chain stdout $f.bigpsl.unmap | sort -k1,1 -k2,2n > $f.bigpsl 
done

###########################
# liftOver maker cdna2genome track
###########################
mkdir big/cdna2genome
mkdir bbi/cdna2genome
cat ../carAur01/big/carAur01.cdna2genome.bigpsl | liftOver -tab -bedPlus=12 stdin $chain stdout big/cdna2genome/carAur03.noM.cdna2genome.bigpsl.unmap | sort -k1,1 -k2,2n > big/cdna2genome/carAur03.noM.cdna2genome.bigpsl
cat ../carAur01/big/carAur01.cdna2genome.ZF.bigpsl | liftOver -tab -bedPlus=12 stdin $chain stdout big/cdna2genome/carAur03.noM.cdna2genome.ZF.bigpsl.unmap | sort -k1,1 -k2,2n > big/cdna2genome/carAur03.noM.cdna2genome.ZF.bigpsl

######################################################
# SNP and DIV from 11549471, 11549472(assembly) and 13778779(WT)
######################################################
mkdir bbi/variant
mkdir big/variant
# copy snv and div
cp /cluster/ifs/users/mullikin/assemblies/goldfish/11549472/fastq/carAur03-2x250/GF72.mpg.???.vcf.gz bbi/variant/
cp /cluster/ifs/users/mullikin/assemblies/goldfish/11549471/fastq/carAur03-2x250/GF71.mpg.???.vcf.gz bbi/variant/
cp /cluster/ifs/users/mullikin/assemblies/goldfish/13778779/fastq/carAur03-2x250/WTgf.mpg.???.vcf.gz bbi/variant/
for f in `ls bbi/variant/*.vcf.gz`
do
    gunzip -c $f | bgzip -c > a; mv a $f;
    tabix -p vcf $f
done
zcat bbi/variant/*.vcf.gz |  awk -F$'\t' -v OFS=$'\t' '$1!~/^#/ {l=length($4)-1; print $1,$2-1,$2+l}' | sort -k1,1 -k2,2n -s 2G | uniq > a1
cat carAur03.sm.fa.fai | awk -F$'\t' -v OFS=$'\t' '{print $1,0,$2}' | sort -k1,1 -k2,2n > a2.bed
bedtools subtract -a a2.bed -b a1.bed | awk -F$'\t' -v l=1000 '$3-$2>=l' > big/variant/invariant.L1000.bed

####################################################
# chain and net
####################################################
mkdir big/chain_net
mkdir bbi/chain_net
for sp in carAur03 danRer10 carp_ncbi
do
mkdir big/chain_net/carAur03.vs.$sp
cd big/chain_net/carAur03.vs.$sp
hgLoadChain -noBin -test carAur03 carAur03_$sp all.target.syn.net.chain
sed 's/\.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > net.bigchain
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > net.biglink
hgLoadChain -noBin -test carAur03 carAur03_$sp all.prenet.chain
sed 's/.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > prenet.bigchain
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > prenet.biglink
cd ../../../
done
sp=danRer10
hgLoadChain -noBin -test carAur03 carAur03_$sp all.prenet.query.chain
sed 's/.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > prenet.query.bigchain
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > prenet.query.biglink

######################################################
# 4-way, 3-way phastcons
######################################################
# {{{
mkdir -p big/phastcons/run2.roast_multiz
for sp in ZF GC CC GF
for a in ZF_GC_CC_GF ZF_GC_CC ZF_GC_GF
do 
do 
f=../chain_net1/multi_run2/$a.out/all.sing.roast_multiz.f.$sp.phastcons.bed ;
if [ -f $f ]; then cp $f big/phastcons/run2.roast_multiz/$a.f.$sp.bed;  fi
f=../chain_net1/multi_run2/$a.out/all.sing.roast_multiz.$sp.phastcons.bed ;
if [ -f $f ]; then cp $f big/phastcons/run2.roast_multiz/$a.$sp.bed;  fi
done
done
#
mkdir -p big/phastcons/run2.net_roast_multic
for sp in ZF GC CC GF
for a in ZF_GC_CC_GF ZF_GC_CC ZF_GC_GF
do 
do 
f=../chain_net1/multi_run2/$a.netSubset.out/all.sing.roast_multic.f.$sp.phastcons.bed ;
if [ -f $f ]; then cp $f big/phastcons/run2.net_roast_multic/$a.f.$sp.bed;  fi
f=../chain_net1/multi_run2/$a.netSubset.out/all.sing.roast_multic.$sp.phastcons.bed ;
if [ -f $f ]; then cp $f big/phastcons/run2.net_roast_multic/$a.$sp.bed;  fi
done
#
done
for sp in ZF GC CC GF
for a in ZF_GC_CC_GF ZF_GC_CC ZF_GC_GF ZF_GF ZF_GC ZF_CC GC_GF CC_GF GC_CC
do 
do 
f=../chain_net1/multi_run2/$a.netSubset.out/all.orig.roast_multic_$sp.$sp.phastcons.bed ;
if [ -f $f ]; then cp $f big/phastcons/run2.net_roast_multic/$a.$sp.bed;  fi
done
done

bedtools intersect -a big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.GF.bed -b big/carAur03.repeatmasker.merged.bed3 -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.GF.notR.bed
cat big/carAur03.noM.exon.bed big/carAur03.noM.ncrna.f.not_ovl.exon.bed |  bedtools intersect -a big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.GF.notR.bed -b stdin -wao | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5 && l<20) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$13;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.GF.notRE.bed

bedtools intersect -a big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.CC.bed -b /home/chenz11/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.repeats.bed -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.CC.notR.bed

bedtools intersect -a big/phastcons/run2.net_roast_multic/CC_GF.CC.bed -b /home/chenz11/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.repeats.bed -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/CC_GF.CC.notR.bed


bedtools intersect -a big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.ZF.bed -b /home/chenz11/data/zebrafish/ensembl85/danRer10.rmsk.bed -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.ZF.notR.bed

bedtools intersect -a big/phastcons/run2.net_roast_multic/ZF_GF.ZF.bed -b /home/chenz11/data/zebrafish/ensembl85/danRer10.rmsk.bed -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/ZF_GF.ZF.notR.bed

for sp in CC GF GC ZF
do
bedtools intersect -a big/phastcons/run2.net_roast_multic/${sp}_GF.GF.bed -b big/carAur03.repeatmasker.merged.bed3 -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/${sp}_GF.GF.notR.bed

cat big/carAur03.noM.exon.bed big/carAur03.noM.ncrna.f.not_ovl.exon.bed |  bedtools intersect -a big/phastcons/run2.net_roast_multic/${sp}_GF.GF.notR.bed -b stdin -wao | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5 && l<20) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$13;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/${sp}_GF.GF.notRE.bed
done

bedtools intersect -a big/phastcons/run2.net_roast_multic/CC_GF.CC.bed -b big/carAur03.repeatmasker.merged.bed3 -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/CC_GF.CC.notR.bed

bedtools intersect -a big/phastcons/run2.net_roast_multic/GC_CC.CC.bed -b big/carAur03.repeatmasker.merged.bed3 -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/GC_CC.CC.notR.bed

bedtools intersect -a big/phastcons/run2.net_roast_multic/ZF_CC.CC.bed -b big/carAur03.repeatmasker.merged.bed3 -wao | sort -k1,1 -k2,2n | awk -v OFS=$'\t' -F$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6} $7!="." {
if (!(chr==$1 && b==$2 && e==$3) ) {
	if (l<(e-b)*0.5) { print chr,b,e,c4,c5,c6; }
	chr=$1; b=$2; e=$3; c4=$4; c5=$5; c6=$6; l=0;
}
l+=$10;
}
END{ print chr,b,e,c4,c5,c6; } ' > big/phastcons/run2.net_roast_multic/ZF_CC.CC.notR.bed
#
# }}}
##################################################

cat big/carAur03.repeatmasker.bed big/carAur03.noM.exon.bed big/carAur03.noM.ncrna.f.not_ovl.bgp | cut -f 1-6 | sort -k1,1 -k2,2n > big/GF.repeat_exon.bed
bedtools merge -i big/GF.repeat_exon.bed > big/GF.repeat_exon.merged.bed 
cat ~/data/danRer10/danRer10.rmsk.bed ~/data/danRer10/danRer10.ensGene.exons.bed | sort -k1,1 -k2,2n > ZF.repeat_exon.bed
cp /home/chenz11/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.repeat_exon*.bed ./big/
cp ~/data/grass_carp/ng/Ctenopharyngodon_idellus_female/C_idella_female.LG.v1.repeat_exon.bed big/GC.repeat_exon.bed
cp ~/data/grass_carp/ng/Ctenopharyngodon_idellus_female/C_idella_female.LG.v1.repeat_exon.merged.bed big/GC.repeat_exon.merged.bed

bedtools subtract -a ZF_GC_GF.f.ZF.bed -b ZF_GC_CC_GF.f.ZF.bed -f 0.1 -A | sort -k1,1 -k2,2n> a.bed
bedtools subtract -a a.bed -b ZF_GC_CC.f.ZF.bed -f 0.1 -A | sort -k1,1 -k2,2n> ZF_GC_GF.f.ZF.no_CC.bed
bedtools subtract -a ZF_GC_GF.f.ZF.no_CC.bed -b ~/data/danRer10/danRer10.ensGene.exons.bed -f 0.1 -A | sort -k1,1 -k2,2n> ZF_GC_GF.f.ZF.no_CC.not_exon.bed

bedtools subtract -a ZF_GC_CC.f.ZF.bed -b ZF_GC_CC_GF.f.ZF.bed -f 0.1 -A | sort -k1,1 -k2,2n> a.bed
bedtools subtract -a a.bed -b ZF_GC_GF.f.ZF.bed -f 0.1 -A | sort -k1,1 -k2,2n> ZF_GC_CC.f.ZF.no_GF.bed
bedtools subtract -a ZF_GC_CC.f.ZF.no_GF.bed -b ~/data/danRer10/danRer10.ensGene.exons.bed -f 0.1 -A | sort -k1,1 -k2,2n> ZF_GC_CC.f.ZF.no_GF.not_exon.bed

bedtools subtract -a ZF_GC_GF.f.GF.bed -b ZF_GC_CC_GF.f.GF.bed -f 0.1 -A | bedtools subtract -a stdin -b GC_CC_GF.f.GF.bed -f 0.1 -A | sort -k1,1 -k2,2n> ZF_GC_GF.f.GF.no_CC.bed
bedtools subtract -a ZF_GC_CC.f.CC.bed -b ZF_GC_CC_GF.f.CC.bed -f 0.1 -A | bedtools subtract -a stdin -b GC_CC_GF.f.CC.bed -f 0.1 -A | sort -k1,1 -k2,2n> ZF_GC_CC.f.CC.no_GF.bed
bedtools subtract -a ZF_GC_CC.f.CC.no_GF.bed -b CC.repeat_exon.bed -A -f 0.1 > ZF_GC_CC.f.CC.no_GF.not_exon.bed 
bedtools subtract -a ZF_GC_GF.f.GF.bed -b ZF_GC_CC_GF.f.GF.bed -f 0.1 -A | bedtools subtract -a stdin -b GC_CC_GF.f.GF.bed -f 0.1 -A | sort -k1,1 -k2,2n> ZF_GC_GF.f.GF.no_CC.bed
bedtools subtract -a ZF_GC_GF.f.GF.no_CC.bed -b GF.repeat_exon.bed -A -f 0.1 > ZF_GC_GF.f.GF.no_CC.not_exon.bed 



exit 0;
