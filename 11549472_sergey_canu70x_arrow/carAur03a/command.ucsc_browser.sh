cwd=`pwd`
datadir=`echo $cwd | sed 's/data_cmd/data/'`

module load ucsc bedtools

cd $datadir

kent=/data/genome/jksrc_v352/kent/src/hg/lib/

asm=carAur03
asm_dir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur03
ucsc_asm_dir=~/data/datashare/fishCarAur1/$asm

asm1=carAur01
asm1_dir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01
asm2=carAur03
asm2_dir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur03

chain=$asm2_dir/carAur01_to_$asm2.liftOver.chain

mkdir -p $asm_dir/bbi
# transform to 2bit
if ! [ -f $asm.chromSizes ]; then faSize -detailed $asm.fa > $asm.chromSizes; fi
if ! [ -f $asm.2bit ]; then faToTwoBit $asm.fa $asm.2bit; fi
if ! [ -f $ucsc_asm_dir/$asm.2bit ]; cp $asm_dir/$asm.2bit $ucsc_asm_dir/
cp $asm.chromSizes $asm_dir/

# add chrM to agp


# GC track
mkdir bbi;
if ! [ -d $ucsc_asm_dir/bbi ]; then mkdir $ucsc_asm_dir/bbi; fi
if ! [ -f $asm_dir/gc.bw ]; then
# hgGcPercent -wigOut -doGaps -file=stdout -win=100 -verbose=0 carAur02 carAur02.2bit | wigEncode stdin gc.wig gc.wib
hgGcPercent -wigOut -doGaps -file=stdout -win=100 -verbose=0 $asm $asm.2bit > bbi/gc.wiggle
wigToBigWig bbi/gc.wiggle $asm.chromSizes bbi/gc.bw
fi
cat $asm_dir/$asm.fa.fai | awk '{print $1"\t0\t"$2}' > $asm_dir/$asm.scf.bed


###########################
# contig stats , add LG
cat $chain | grep '^chain' | cut -d" " -f 3,8-12 | sort -k1,1 | sed 's/ /\t/g' > a
join -t$'\t' -j 1 ../carAur01/carAur01.ctgInfo2  a > b1
join -t$'\t' -j 1 -v 1 ../carAur01/carAur01.ctgInfo2  a | awk '{print $0"\t"$1"\t"$2"\t+\t0\t"$2}'> b2
cat b1 b2 | sort -k1,1 > carAur01.ctgInfo3
echo "`head -n 1 ../../asm.contigs.layout.tigInfo`"$'\t'"relative_depth"$'\t'"LG"$'\t'"LG_size"$'\t'"strand"$'\t'"begin"$'\t'"end" > carAur01.ctgInfo3.header.txt

>carAur01.ctgInfo3.stat
cat carAur01.ctgInfo3 | awk -F$'\t' -v OFS=$'\t' -v n0=0 -v n1=0 -v n2=0 -v m0=0 -v m1=0 -v m2=0 '$11~/^LG/ {if ($10<0.6) {n0++; m0+=$2;} else if ($10<1.6) {n1++;m1+=$2;} else {n2++;m2+=$2;}} END{print "Coverage_of_LG_contig","Half","One","Multiple"; print "N_CTG",n0,n1,n2; print "total_bp",m0,m1,m2}' > carAur01.ctgInfo3.stat

echo >>  carAur01.ctgInfo3.stat
cat carAur01.ctgInfo3 | awk -F$'\t' -v OFS=$'\t' -v n0=0 -v n1=0 -v m0=0 -v m1=0 '$11~/^LG/ {if ($7=="no") {n0++; m0+=$2;} else {n1++;m1+=$2;} } END{print "Type","Norm","Repeat"; print "N_CTG",n0,n1; print "total_bp",m0,m1}' >> carAur01.ctgInfo3.stat

echo >>  carAur01.ctgInfo3.stat
echo "Half_coverage LG counts\n" >>  carAur01.ctgInfo3.stat;
echo "LG"$'\t'"LG_size"$'\t'"ctg_N"$'\t'"ctg_bp" >>  carAur01.ctgInfo3.stat; 
cat carAur01.ctgInfo3 | perl -e 'while(<>) {chomp; @t=split "\t"; if ($t[10]=~m/^LG/ && $t[9]<0.6) {
$count{$t[10]}[0]=$t[11]; $count{$t[10]}[1]++; $count{$t[10]}[2]+=$t[1];}}
foreach my $chr (sort keys(%count)) {print "$chr\t$count{$chr}[0]\t$count{$chr}[1]\t$count{$chr}[2]\n";}'  >>  carAur01.ctgInfo3.stat; 

###########################

###########################
# repeatmasker
###########################
liftOver -multiple  -bedPlus=12 -tab $asm1_dir/$asm1.repeatmasker.f.bgp $chain stdout $asm2_dir/$asm2.repeatmasker.bgp.unmap | sort -k1,1 -k2,2n > $asm2_dir/$asm2.repeatmasker.bgp
bedToBigBed -type=bed12+8 -tab -as=$kent/bigGenePred.as -extraIndex=name,name2,geneName,geneName2 $asm2_dir/$asm2.repeatmasker.bgp $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$asm2.repeatmasker.bb
cat $asm2_dir/$asm2.repeatmasker.bgp | cut -f 4,13 > $asm2_dir/$asm2.repeatmasker.attr.Name
ixIxx $asm2_dir/$asm2.repeatmasker.attr.Name $asm2_dir/bbi/$asm2.repeatmasker.ix $asm2_dir/bbi/$asm2.repeatmasker.ixx
cat $asm2_dir/$asm2.repeatmasker.bgp | awk -v OFS=$'\t' -F$'\t' '{print $1,$2,$3,$13,$5,$6}' > $asm2_dir/$asm2.repeatmasker.bed
# ----------------------------
# repeat bed
liftOver -multiple  -bedPlus=12 -tab $asm1_dir/$asm1.repeat.bed $chain stdout $asm2_dir/$asm2.repeat.bed.unmap | bedtools sort -i > $asm2_dir/$asm2.repeat.bed
bedtools merge -i carAur03.repeat.bed -c 4,5 -o distinct,count > carAur03.repeat.bedtools_merge.bed
bedtools maskfasta -soft -fi carAur03.fa -fo carAur03.sm.fa -bed carAur03.repeat.bed
faToTwoBit $asm2.sm.fa $asm2.2bit
samtools faidx carAur03.sm.fa
~/my_program3/src/assembler/NL50.pl carAur03.sm.fa.fai > carAur03.sm.fa.NL50
cat carAur03.sm.fa.fai | cut -f 1,2 | grep arrow | sort -k1,1 > carAur03.unplace.names

###########################
# GC content
###########################
#gccontent carAur03.sm.fa > carAur03.gccontent.txt
#bedtools nuc -fi carAur03.sm.fa -bed carAur03.scf.bed > carAur03.bedtools_nuc_out.txt
~/my_program3/src/utility/czl_fasta_stat -i carAur03.sm.fa -o carAur03.fasta_stat_out > carAur03.fasta_stat_out.stdout

#######################################################
# gene track
#######################################################
liftOver -multiple -genePred $asm1_dir/$asm1.gene.annot.gp $chain $asm2_dir/$asm2.gene.annot.gp $asm2_dir/$asm2.gene.annot.gp.unmap
genePredToBigGenePred $asm2.gene.annot.gp stdout | sort -k1,1 -k2,2n > $asm2_dir/$asm2.gene.annot.bgp
bedToBigBed -type=bed12+8 -tab -as=$kent/bigGenePred.as -extraIndex=name,name2,geneName,geneName2 $asm2_dir/$asm2.gene.annot.bgp $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$asm2.gene.annot.bb
cat $asm2_dir/$asm2.gene.annot.bgp | cut -f 4,13 > $asm2_dir/$asm2.gene.annot.attr.Name
ixIxx $asm2_dir/$asm2.gene.annot.attr.Name $asm2_dir/bbi/$asm2.gene.annot.ix $asm2_dir/bbi/$asm2.gene.annot.ixx
cut -f 1-6 $asm2_dir/$asm2.gene.annot.bgp > $asm2_dir/$asm2.gene.annot.bed
liftOver -multiple $asm1_dir/$asm1.gene.annot.exon.bed $chain $asm2_dir/$asm2.gene.annot.exon.bed $asm2_dir/$asm2.gene.annot.exon.unmap
##
# produce carAur03.CDS.for_bigpsl
cat carAur03.gene.annot.bgp | perl -ne '
if (m/^#/ || m/^\s*$/) { next; }
chomp;
my @t=split "\t",$_;  my $b; my $e;
$cds_start=$t[6];
$cds_end=$t[7];
if ($t[9]==1) {
    $b = $t[1]+$t[11];
    $e = $t[1]+$t[11]+$t[10];
    if ($b<$cds_start) { $b = $cds_start; }
    if ($e>$cds_end) { $e = $cds_end; }

    print $t[0] , "\t" , $b+1 , ".." , $e , "\n";
} else {
    my @sz = split ",", $t[10];
    my @b = split ",", $t[11];
    my @aa = ();
    for (my $i=0; $i<=$#sz; $i++) {
        $b = $t[1]+$b[$i];
        $e = $b+$sz[$i];
        if ($e <= $cds_start) {next;}
        if ($b >= $cds_end) {next;}
        if ($b<$cds_start) { $b = $cds_start; }
        if ($e>$cds_end) { $e = $cds_end; }
        push @aa, $b+1 . ".." . $e; 
    }
	if ($#aa==0) { print "$t[0]\t$aa[0]\n"; }
	else { print "$t[0]\tjoin(", join(",", @aa), ")", "\n"; }
}' > carAur03.CDS.for_bigpsl

#######################################################
# est2genome 
#######################################################
liftOver -multiple -pslT $asm1_dir/$asm1.est2genome.psl $chain $asm2_dir/$asm2.est2genome.psl $asm2_dir/$asm2.est2genome.psl.unmap
pslToBigPsl $asm2_dir/$asm2.est2genome.psl stdout | sort -k1,1 -k2,2n > $asm2.est2genome.bigpsl
bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=$kent/bigPsl.as "$asm2_dir/$asm2.est2genome.bigpsl" "$asm2_dir/$asm2.chromSizes" "$asm2_dir/bbi/$asm2.est2genome.bb"

#######################################################
# cdna2genome 
#######################################################
liftOver -bedPlus=12 -tab $asm1_dir/$asm1.cdna2genome.bigpsl $chain stdout $asm2_dir/$asm2.cdna2genome.bigpsl.unmap  | sort -k1,1 -k2,2n > $asm2_dir/$asm2.cdna2genome.bigpsl
bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=$kent/bigPsl.as "$asm2_dir/$asm2.cdna2genome.bigpsl" "$asm2_dir/$asm2.chromSizes" "$asm2_dir/bbi/$asm2.cdna2genome.bb"
cp $asm2_dir/bbi/$asm2.cdna2genome.bigpsl.bb ~/data/datashare/fishCarAur1/carAur03/bbi/

#######################################################
# cdna2genome track
echo "##gff-version 3" >$asm1.gene.cdna2genome.gff
zcat ../maker4a/carAur01.gene.cdna2genome.gff.gz | perl -ne 'if (m/^#/) { print; next; }
    chomp;
    my @tab = split "\t";
    my $s = ($tab[6] eq "+") ? "F" : "R";
    my %default=(ID=>0, Name=>0, Alias=>0, Parent=>0, Target=>0, Gap=>0, Derives_from=>0, Note=>0, Xref=>0, Ontology_term=>0);
    my $attrs="";
    if ($tab[2] eq "expressed_sequence_match") {$tab[2]="cDNA_match"; }
    foreach $attr (split /\s*;\s*/, $tab[8]) {
        my ($u,$v) = split "=", $attr;
        if (exists $default{$u}) {
        } else {
            $attr="u$attr";
        }
        if ($attrs ne "") { $attrs.=";"; }
        $attrs .= $attr;
    }
    $tab[8] = $attrs;
    print join("\t",@tab), "\n";
    ' >> $asm1.gene.cdna2genome.gff
perl ~/my_program3/src/utility/czl_gff_match_part_to_match.pl -i $asm1.gene.cdna2genome.gff -o $asm1.gene.cdna2genome.1.gff
gff3ToPsl $genome1.chromSizes ~/data/ensembl85/ens85_and_carp_gcarp.rna.short_name.chromSizes $asm1.gene.cdna2genome.1.gff $asm1.gene.cdna2genome.psl
pslSwap $asm1.gene.cdna2genome.psl $asm1.gene.cdna2genome.psl.tmp
mv $asm1.gene.cdna2genome.psl.tmp $asm1.gene.cdna2genome.psl
pslPosTarget $asm1.gene.cdna2genome.psl $asm1.gene.cdna2genome.psl.tmp 
mv $asm1.gene.cdna2genome.psl.tmp $asm1.gene.cdna2genome.psl
# add gene name
cat $asm1.gene.cdna2genome.psl | awk '{print $10"\t"$0}' | sort -k1,1 > $asm1.gene.cdna2genome.psl.tmp
sort -k1,1 ~/data/ensembl85/ens85.tid_cdna_ttype_gid_gname_gtype.danRer10 > $asm1.gene.cdna2genome.psl.tmp2
join -t $'\t' -j 1 $asm1.gene.cdna2genome.psl.tmp $asm1.gene.cdna2genome.psl.tmp2 | awk '{a=$2; for (i=3;i<=10; i++) {a=a"\t"$i;} a=a"\t"$11"__"$26; for (i=12;i<=22;i++) {a=a"\t"$i;} print a}' > $asm1.gene.cdna2genome.psl.tmp3 
join -t $'\t' -j 1 -v 1 $asm1.gene.cdna2genome.psl.tmp $asm1.gene.cdna2genome.psl.tmp2 | cut -f 2-22 >> $asm1.gene.cdna2genome.psl.tmp3 
cat $asm1.gene.cdna2genome.psl.tmp3 > $asm1.gene.cdna2genome.psl
rm $asm1.gene.cdna2genome.psl.tmp*
#
pslToBigPsl $asm1.gene.cdna2genome.psl stdout | sort -k1,1 -k2,2n > $asm1.gene.cdna2genome.bigpsl
bedToBigBed -type=bed12+13 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigPsl.as $asm1.gene.cdna2genome.bigpsl $genome1.chromSizes $asm1.gene.cdna2genome.bb
#
liftOver -pslT $asm1.gene.cdna2genome.psl $chain $asm2.gene.cdna2genome.psl $asm2.gene.cdna2genome.psl.unmap
pslToBigPsl $asm2.gene.cdna2genome.psl stdout | sort -k1,1 -k2,2n > $asm2.gene.cdna2genome.bigpsl
bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigPsl.as $asm2.gene.cdna2genome.bigpsl $genome2.chromSizes $asm2.gene.cdna2genome.bb

###########################
# tRNA
###########################
liftOver -multiple -genePred $asm1_dir/$asm1.tRNA.gp $chain $asm2_dir/$asm2.tRNA.gp $asm2_dir/$asm2.tRNA.gp.unmap
genePredToBigGenePred $asm2.tRNA.gp stdout | sort -k1,1 -k2,2n > $asm2_dir/$asm2.tRNA.bgp
bedToBigBed -type=bed12+8 -tab -as=$kent/bigGenePred.as -extraIndex=name,name2,geneName,geneName2 $asm2_dir/$asm2.tRNA.bgp $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$asm2.tRNA.bb
cat $asm2_dir/$asm2.tRNA.bgp | cut -f 4,13 > $asm2_dir/$asm2.tRNA.attr.Name
ixIxx $asm2_dir/$asm2.tRNA.attr.Name $asm2_dir/bbi/$asm2.tRNA.ix $asm2_dir/bbi/$asm2.tRNA.ixx


###########################
# miRNA
###########################
fn=$asm2.miRNA
liftOver -bedPlus=12 -tab $asm1_dir/$asm1.miRNA.bigpsl $chain stdout $asm2_dir/$asm2.miRNA.bigpsl.unmap  | sort -k1,1 -k2,2n > $asm2_dir/$asm2.miRNA.bigpsl
bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigPsl.as $asm2_dir/$asm2.miRNA.bigpsl $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$asm2.miRNA.bigpsl.bb
cp $asm2_dir/bbi/$asm2.miRNA.bigpsl.bb ~/data/datashare/fishCarAur1/carAur03/bbi/


###########################
# Rfam
###########################
fn=$asm2.rfam
liftOver -bedPlus=12 -tab $asm1_dir/$asm1.rfam.bgp $chain stdout $asm2_dir/$asm2.rfam.bgp.unmap | sort -k1,1 -k2,2n > $asm2_dir/$asm2.rfam.bgp
bedToBigBed -extraIndex=name,name2 -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as $asm2_dir/$asm2.rfam.bgp $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$asm2.rfam.bgp.bb
cp $asm2_dir/bbi/$asm2.rfam.bgp.bb ~/data/datashare/fishCarAur1/carAur02/bbi/


############################################
# phastcons and conserved region
############################################
# region
fn=$asm2.roast.mz.phastcon.GC_GF_CC
liftOver -tab ../chain_net1/multiz4/gcarpF.carpNG.carAur/roast.mz.carAur.bed $chain bbi/$fn.bed bbi/$fn.bed.unmap
bedToBigBed bbi/$fn.bed $asm2.chromSizes bbi/$fn.bb
fn=$asm2.tba.mz.phastcon.GC_GF_CC
liftOver -tab ../chain_net1/multiz4/gcarpF.carpNG.carAur/tba.mz.carAur.bed $chain bbi/$fn.bed bbi/$fn.bed.unmap
bedToBigBed bbi/$fn.bed $asm2.chromSizes bbi/$fn.bb

fn=$asm2.ZF_GC_GF_CC.roast.mz.phastcons
liftOver $asm1_dir/$asm1.ZF_GC_GF_CC.roast.mz.phastcons.bed $chain stdout $fn.bed.unmap | sort -k1,1 -k2,2n > $asm2_dir/$fn.bed
bedToBigBed $asm2_dir/$fn.bed $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$fn.bb
fn=$asm2.ZF_GC_GF_CC.tba.mz.phastcons
liftOver $asm1_dir/$asm1.ZF_GC_GF_CC.tba.mz.phastcons.bed $chain stdout $fn.bed.unmap | sort -k1,1 -k2,2n > $asm2_dir/$fn.bed
bedToBigBed $asm2_dir/$fn.bed $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$fn.bb

# phastcons peaks 
fn=$asm2.ZF_GC_GF_CC.roast.mz.phastcons
bigWigToBedGraph $asm1_dir/bbi/$asm1.ZF_GC_GF_CC.roast.mz.phastcons.bw stdout | liftOver stdin $chain stdout $asm2_dir/$fn.bedgraph.unmap | sort -k1,1 -k2,2g > $asm2_dir/$fn.bg
bedGraphToBigWig $asm2_dir/$fn.bg $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$fn.bw
fn=$asm2.ZF_GC_GF_CC.tba.mz.phastcons
bigWigToBedGraph $asm1_dir/bbi/$asm1.ZF_GC_GF_CC.tba.mz.phastcons.bw stdout | liftOver stdin $chain stdout $asm2_dir/$fn.bedgraph.unmap | sort -k1,1 -k2,2g > $asm2_dir/$fn.bg
bedGraphToBigWig $asm2_dir/$fn.bg $asm2_dir/$asm2.chromSizes $asm2_dir/bbi/$fn.bw


############################################
# chains and nets
############################################
zcat $asm2_dir/../chain_net1/carAur.vs.carAur/all.prenet.fself.chain.gz | chainToPsl stdin $asm2_dir/$asm2.chromSizes $asm2_dir/$asm2.chromSizes $genome1fa $genome1fa STDOUT | gzip -c > ../chain_net1/carAur.vs.carAur/all.prenet.chain.psl.gz
liftOver -pslT ../chain_net1/carAur.vs.carAur/all.prenet.chain.psl $chain $asm2.self.chain.psl $asm2.self.chain.psl.unmap
# pslToChain $asm2.self.chain.psl $asm2.self.chain
pslSwap $asm2.self.chain.psl $asm2.self.chain.swap.psl
liftOver -pslT $asm2.self.chain.swap.psl $chain $asm2.self.chain.psl $asm2.self.chain.psl.Q.unmap
pslToChain $asm2.self.chain.psl $asm2.self.chain
hgLoadChain -noBin -test $asm2 bigChain self $asm2.self.chain
sed 's/.000000//' $asm2.self.chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $asm2.self.bigchain
bedToBigBed -type=bed6+6 -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigChain.as -tab $asm2.self.bigchain $asm2.chromSizes bbi/$asm2.self.chain.bb

if ! [ -f ../chain_net1/carAur.vs.carAur/all.target.net.chain.psl.gz ]
then
    chainToPsl ../chain_net1/carAur.vs.carAur/all.target.net.chain $genome1.fa.fai $genome1.fa.fai $genome1fa $genome1fa STDOUT | gzip -c > ../chain_net1/carAur.vs.carAur/all.target.net.chain.psl.gz
fi
liftOver -pslT ../chain_net1/carAur.vs.carAur/all.prenet.psl $chain $asm2.self.chain.psl $asm2.self.chain.psl.unmap
pslToChain $asm2.self.chain.psl $asm2.self.chain
hgLoadChain -noBin -test $asm2 bigChain self $asm2.self.chain
sed 's/.000000//' $asm2.self.chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $asm2.self.bigchain
bedToBigBed -type=bed6+6 -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigChain.as -tab $asm2.self.bigchain $asm2.chromSizes bbi/$asm2.self.chain.bb

perl ~/my_program3/src/utility/czl_liftover.pl -i ../chain_net1/carAur.vs.carAur/all.target.net.chain -L $chain -LT -o - | chainSwap /dev/stdin /dev/stdout | perl ~/my_program3/src/utility/czl_liftover.pl -i - -L $chain -LT -o  $asm2.target.net.chain.gz
zcat $asm2.self.net.chain.gz | hgLoadChain -noBin -test $asm2 self /dev/stdin 
mv chain.tab $asm2.self.net.chain.tab
mv link.tab $asm2.self.net.link.tab
sed 's/.000000//' $asm2.self.net.chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $asm2.self.net.bigchain
bedToBigBed -type=bed6+6 -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigChain.as -tab $asm2.self.net.bigchain $asm2.chromSizes bbi/$asm2.self.net.chain.bb
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' $asm2.self.net.link.tab | sort -k1,1 -k2,2n >$asm2.self.net.biglink
bedToBigBed -type=bed4+1 -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigLink.as -tab $asm2.self.net.biglink $asm2.chrom.sizes bbi/$asm2.self.net.chain.link.bb

zcat ../chain_net1/carAur.vs.carAur/all.prenet.chain | chainSplit -lump=200 ../chain_net1/carAur.vs.carAur/all.prenet.chain_split 

~/my_program3/src/utility chain_to_bb.sh carAur03.danRer10.target.net.chain.gz carAur03.fa.fai bbi/carAur03.danRer10.target.net
~/my_program3/src/utility chain_to_bb.sh carAur03.danRer10.prenet.chain.gz carAur03.fa.fai bbi/carAur03.danRer10.prenet

liftOver -pslT ../chain_net1/carAur.vs.danRer10/all.prenet.psl $chain $asm2.vs.danRer10.net.chain.psl $asm2.self.net.chain.psl.unmap
liftOver -pslT ../chain_net1/carAur.vs.danRer10/all.target.net.psl $chain $asm2.vs.danRer10.net.chain.psl $asm2.self.net.chain.psl.unmap
pslPosTarget $asm2.vs.danRer10.net.chain.psl $asm2.vs.danRer10.net.chain.psl.tmp; mv  $asm2.vs.danRer10.net.chain.psl.tmp  $asm2.vs.danRer10.net.chain.psl

# chain
# This is WRONG
# liftOver -bedPlus=12 -tab $asm1_dir/$asm1.vs.danRer10.net.bigchain $chain stdout $asm2_dir/$asm2.vs.danRer10.net.bigchain.unmap | sort -k1,1 -k2,2n > $asm2_dir/$asm2.vs.danRer10.net.bigchain

cat carAur03.gene.annot.bgp | perl -ne '@t=split "\t",$_,-1; my $g=$t[3]; $g=~s/-mRNA-[0-9]+$//; $l=0; foreach my $l1 (split(",", $t[10])) {$l+=$l1;} print "$t[0]\t$t[1]\t$t[2]\tcarAur|$g|$t[3]|$t[3]|$t[12]|$l|$l\t$t[4]\t$t[5]\n";' > carAur03.gene.annot.6.bed

# copy to ucsc folder
files=$asm.2bit bbi/gc.bw
for f in $files
do
    cp $asm_dir/$f $ucsc_asm_dir/
    done

cd $cwd


