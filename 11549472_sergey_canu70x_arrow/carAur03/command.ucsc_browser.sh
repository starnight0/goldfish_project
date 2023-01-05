cwd=`pwd`
datadir=`echo $cwd | sed 's/data_cmd/data/'`

module load ucsc bedtools

cd $datadir

kent=/data/genome/jksrc_v352/kent/src/hg/lib/

asm=carAur03
asm_dir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur03
ucsc_dir=~/data/datashare/fishCarAur1/carAur03
mkdir $ucsc_dir/bbi
mkdir $ucsc_dir/bbi/est2genome

# assembly
bedToBigBed big/carAur03.assembly.bed carAur03.sm.fa.fai bbi/carAur03.assembly.bb;
cp bbi/carAur03.assembly.bb ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.assembly.bb;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.assembly.bb;

# repeatmasker
bedToBigBed -extraIndex=name big/carAur03.repeatmasker.bed carAur03.sm.fa.fai bbi/carAur03.repeatmasker.bb;
cp bbi/carAur03.repeatmasker.bb ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.repeatmasker.bb;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.repeatmasker.bb;

# gene
cat big/carAur03.noM.gene.bgp | cut -f 4,19  > a
#cat big/carAur03.noM.gene.bgp | cut -f 4,18 >> a
#cat big/carAur03.noM.gene.bgp | cut -f 4,19 >> a
#cat big/carAur03.noM.gene.bgp | awk -F$'\t' -v OFS=$'\t' '{$4=$4" "$19; print $0}' > b;
#bedToBigBed -type=bed12+8 -tab -as=$kent/bigGenePred.as -extraIndex=name,name2,geneName,geneName2 b carAur03.sm.fa.fai bbi/carAur03.gene.bb
bedToBigBed -type=bed12+8 -tab -as=$kent/bigGenePred.as -extraIndex=name,name2,geneName,geneName2 big/carAur03.noM.gene.bgp carAur03.sm.fa.fai bbi/carAur03.gene.bb
ixIxx a bbi/carAur03.gene.ix bbi/carAur03.gene.ixx
cp bbi/carAur03.gene.bb ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.gene.bb ;
cp bbi/carAur03.gene.ix ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.gene.ix ;
cp bbi/carAur03.gene.ixx ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.gene.ixx ;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.gene.*;

bedToBigBed big/carAur03.noM.gene.masked.bed carAur03.sm.fa.fai bbi/carAur03.gene.masked.bb
cp bbi/carAur03.gene.masked.bb ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.gene.masked.bb ;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.gene.masked.bb;

# ncRNA
#cat big/carAur03.noM.ncrna.f.bgp | cut -f 4,18 > a
#cat big/carAur03.noM.ncrna.f.bgp | cut -f 4,19 >> a
bedToBigBed -type=bed12+8 -tab -as=$kent/bigGenePred.as -extraIndex=name,name2,geneName big/carAur03.noM.ncrna.f.not_ovl.bgp carAur03.sm.fa.fai bbi/carAur03.noM.ncrna.f.bb
#ixIxx a bbi/carAur03.noM.ncrna.f.ix bbi/carAur03.noM.ncrna.f.ixx
cp bbi/carAur03.noM.ncrna.f.bb ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.noM.ncrna.f.bb ;
#cp bbi/carAur03.noM.ncrna.f.ix ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.noM.ncrna.f.ix ;
#cp bbi/carAur03.noM.ncrna.f.ixx ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.noM.ncrna.f.ixx ;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.noM.ncrna.*;

# tRNA
bedToBigBed -type=bed12+8 -tab -as=$kent/bigGenePred.as -extraIndex=name,name2,geneName,geneName2 big/carAur03.tRNA.bgp carAur03.sm.fa.fai bbi/carAur03.tRNA.bgp.bb
cp bbi/carAur03.tRNA.bgp.bb ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.tRNA.bgp.bb ;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/bbi/carAur01.tRNA.bgp.bb

# miRNA
f=big/carAur03.noM.miRNA
f1=big/carAur01.noM.miRNA
b=bbi/carAur03.noM.miRNA
b1=bbi/carAur01.noM.miRNA
bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as $f.I95.bigpsl carAur03.sm.fa.fai $b.I95.bigpsl.bb
cp $b.I95.bigpsl.bb ~/data/datashare/fishCarAur1/carAur01/$b1.I95.bigpsl.bb;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/$b1.I95.bigpsl.bb;
#bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as $f.bigpsl carAur03.sm.fa.fai $b.bigpsl.bb
#cp $b.bigpsl.bb ~/data/datashare/fishCarAur1/carAur01/$b1.bigpsl.bb;
#setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/$b1.bigpsl.bb;

# rfam
f=big/carAur03.noM.rfam
f1=big/carAur01.noM.rfam
b=bbi/carAur03.noM.rfam
b1=bbi/carAur01.noM.rfam
bedToBigBed -type=bed12+8 -extraIndex=name,name2,geneName,geneName2 -tab -as=$kent/bigGenePred.as $f.bgp carAur03.sm.fa.fai $b.bgp.bb
cp $b.bgp.bb ~/data/datashare/fishCarAur1/carAur01/$b1.bgp.bb;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/$b1.bgp.bb;

# Interpro
f=big/carAur03.noM.ips
f1=`echo $f | sed 's/carAur03/carAur01/'`
b=bbi/carAur03.noM.ips
b1=`echo $b | sed 's/carAur03/carAur01/'`
bedToBigBed -type=bed12+8 -extraIndex=name,name2,geneName,geneName2 -tab -as=$kent/bigGenePred.as $f.bgp carAur03.sm.fa.fai $b.bgp.bb
cp $b.bgp.bb ~/data/datashare/fishCarAur1/carAur01/$b1.bgp.bb;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/$b1.bgp.bb;
# GO
f=big/carAur03.noM.GO
b=bbi/carAur03.noM.GO
f1=`echo $f | sed 's/carAur03/carAur01/'`
b1=`echo $b | sed 's/carAur03/carAur01/'`
bedToBigBed -type=bed12+8 -extraIndex=name,name2,geneName,geneName2 -tab -as=$kent/bigGenePred.as $f.bgp carAur03.sm.fa.fai $b.bgp.bb
cp $b.bgp.bb ~/data/datashare/fishCarAur1/carAur01/$b1.bgp.bb;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/$b1.bgp.bb;


# est2genome
for sm in `cat ../Trinity_GG_run01/sample | head -n 1`
do
    f=big/est2genome/$sm.trinity_GG
    b=bbi/est2genome/$sm.trinity_GG
    f1=`echo $f | sed 's/carAur03/carAur01/'`
    b1=`echo $b | sed 's/carAur03/carAur01/'`
    chain=carAur01_to_carAur03.liftOver.chain
    bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as $f.bigpsl carAur03.sm.fa.fai $b.bigpsl.bb
    cp $b.bigpsl.bb ~/data/datashare/fishCarAur1/carAur01/$b1.bigpsl.bb
    setfacl -m user:webcpu:r--  ~/data/datashare/fishCarAur1/carAur01/$b1.bigpsl.bb
done

# cdna2genome
mkdir ~/data/datashare/fishCarAur1/carAur03/bbi/cdna2genome
f=big/cdna2genome/carAur03.noM.cdna2genome.ZF
b=bbi/cdna2genome/carAur03.noM.cdna2genome.ZF
f1=`echo $f | sed 's/carAur03/carAur01/'`
b1=`echo $b | sed 's/carAur03/carAur01/'`
bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as $f.bigpsl carAur03.sm.fa.fai $b.bigpsl.bb
cp $b.bigpsl.bb ~/data/datashare/fishCarAur1/carAur01/$b1.bigpsl.bb
setfacl -m user:webcpu:r--  ~/data/datashare/fishCarAur1/carAur01/$b1.bigpsl.bb
#
f=big/cdna2genome/carAur03.noM.cdna2genome
b=bbi/cdna2genome/carAur03.noM.cdna2genome
f1=`echo $f | sed 's/carAur03/carAur01/'`
b1=`echo $b | sed 's/carAur03/carAur01/'`
bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as $f.bigpsl carAur03.sm.fa.fai $b.bigpsl.bb
cp $b.bigpsl.bb ~/data/datashare/fishCarAur1/carAur01/$b1.bigpsl.bb
setfacl -m user:webcpu:r--  ~/data/datashare/fishCarAur1/carAur01/$b1.bigpsl.bb

# variants
cp -r bbi/variant ~/data/datashare/fishCarAur1/carAur03/bbi/
bedToBigBed big/variant/invariant.L1000.bed carAur03.sm.fa.fai bbi/variant/invariant.L1000.bb
cp bbi/variant/invariant.L1000.bb ~/data/datashare/fishCarAur1/carAur03/bbi/variant/
setfacl -m user:webcpu:r-x  ~/data/datashare/fishCarAur1/carAur03/bbi/variant
setfacl -m user:webcpu:r--  ~/data/datashare/fishCarAur1/carAur03/bbi/variant/*
chmod 0640 ~/data/datashare/fishCarAur1/carAur03/bbi/variant/*

###############################################
# chain and net
###############################################
mkdir  ~/data/datashare/fishCarAur1/carAur03/bbi/chain_net
setfacl -m user:webcpu:r-x  ~/data/datashare/fishCarAur1/carAur03/bbi/chain_net
for sp in carAur03 danRer10 carp_ncbi
do
    bedToBigBed -type=bed6+6 -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigChain.as -tab big/chain_net/carAur03.vs.$sp/net.bigchain carAur03.sm.fa.fai bbi/chain_net/carAur03_$sp.net.bigchain.bb
    bedToBigBed -type=bed4+1 -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigLink.as -tab big/chain_net/carAur03.vs.$sp/net.biglink carAur03.sm.fa.fai bbi/chain_net/carAur03_$sp.net.biglink.bb

    bedToBigBed -type=bed6+6 -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigChain.as -tab big/chain_net/carAur03.vs.$sp/prenet.bigchain carAur03.sm.fa.fai bbi/chain_net/carAur03_$sp.prenet.bigchain.bb
    bedToBigBed -type=bed4+1 -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigLink.as -tab big/chain_net/carAur03.vs.$sp/prenet.biglink carAur03.sm.fa.fai bbi/chain_net/carAur03_$sp.prenet.biglink.bb
done
cp bbi/chain_net/*  ~/data/datashare/fishCarAur1/carAur03/bbi/chain_net/
setfacl -m user:webcpu:r-x  ~/data/datashare/fishCarAur1/carAur03/bbi/chain_net/*
chmod 0644 ~/data/datashare/fishCarAur1/carAur03/bbi/chain_net/*

###############################################
# phastcons and conserved region
###############################################
# {{{
mkdir -p $ucsc_dir/bbi/phastcons/run2.roast_multiz
mkdir -p $ucsc_dir/bbi/phastcons/run2.net_roast_multiz
mkdir -p $ucsc_dir/bbi/phastcons/run2.net_roast_multic
setfacl -m user:webcpu:r-x $ucsc_dir/bbi/phastcons/run2.roast_multiz;
setfacl -m user:webcpu:r-x $ucsc_dir/bbi/phastcons/run2.net_roast_multiz;
setfacl -m user:webcpu:r-x $ucsc_dir/bbi/phastcons/run2.net_roast_multic;
cp big/phastcons/run2.roast_multiz/ZF_GC_CC_GF.GF.bb $ucsc_dir/bbi/phastcons/run2.roast_multiz/;
cp big/phastcons/run2.roast_multiz/ZF_GC_CC_GF.GF.bw $ucsc_dir/bbi/phastcons/run2.roast_multiz/;
cp big/phastcons/run2.roast_multiz/ZF_GC_CC_GF.f.GF.bb $ucsc_dir/bbi/phastcons/run2.roast_multiz/;
cp big/phastcons/run2.roast_multiz/ZF_GC_CC_GF.f.GF.bw $ucsc_dir/bbi/phastcons/run2.roast_multiz/;
cp big/phastcons/run2.roast_multiz/ZF_GC_GF.f.GF.bb $ucsc_dir/bbi/phastcons/run2.roast_multiz/;
cp big/phastcons/run2.roast_multiz/ZF_GC_GF.f.GF.bw $ucsc_dir/bbi/phastcons/run2.roast_multiz/;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur03/bbi/phastcons/run2.roast_multiz/*;
# net_roast_multiz
cp big/phastcons/run2.net_roast_multiz/ZF_GC_CC_GF.GF.bb $ucsc_dir/bbi/phastcons/run2.net_roast_multiz/;
cp big/phastcons/run2.net_roast_multiz/ZF_GC_CC_GF.GF.bw $ucsc_dir/bbi/phastcons/run2.net_roast_multiz/;
cp big/phastcons/run2.net_roast_multiz/ZF_GC_CC_GF.f.GF.bb $ucsc_dir/bbi/phastcons/run2.net_roast_multiz/;
cp big/phastcons/run2.net_roast_multiz/ZF_GC_CC_GF.f.GF.bw $ucsc_dir/bbi/phastcons/run2.net_roast_multiz/;
cp big/phastcons/run2.net_roast_multiz/ZF_GC_GF.GF.bb $ucsc_dir/bbi/phastcons/run2.net_roast_multiz/;
cp big/phastcons/run2.net_roast_multiz/ZF_GC_GF.GF.bw $ucsc_dir/bbi/phastcons/run2.net_roast_multiz/;
cp big/phastcons/run2.net_roast_multiz/ZF_GC_GF.f.GF.bb $ucsc_dir/bbi/phastcons/run2.net_roast_multiz/;
cp big/phastcons/run2.net_roast_multiz/ZF_GC_GF.f.GF.bw $ucsc_dir/bbi/phastcons/run2.net_roast_multiz/;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur03/bbi/phastcons/run2.net_roast_multiz/*;

# net_roast_multic
cp big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.GF.bb $ucsc_dir/bbi/phastcons/run2.net_roast_multic/;
cp big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.GF.bw $ucsc_dir/bbi/phastcons/run2.net_roast_multic/;
cp big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.f.GF.bb $ucsc_dir/bbi/phastcons/run2.net_roast_multic/;
cp big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.f.GF.bw $ucsc_dir/bbi/phastcons/run2.net_roast_multic/;
cp big/phastcons/run2.net_roast_multic/ZF_GC_GF.GF.bb $ucsc_dir/bbi/phastcons/run2.net_roast_multic/;
cp big/phastcons/run2.net_roast_multic/ZF_GC_GF.GF.bw $ucsc_dir/bbi/phastcons/run2.net_roast_multic/;
cp big/phastcons/run2.net_roast_multic/ZF_GC_GF.f.GF.bb $ucsc_dir/bbi/phastcons/run2.net_roast_multic/;
cp big/phastcons/run2.net_roast_multic/ZF_GC_GF.f.GF.bw $ucsc_dir/bbi/phastcons/run2.net_roast_multic/;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur03/bbi/phastcons/run2.net_roast_multic/*;
#}}}
