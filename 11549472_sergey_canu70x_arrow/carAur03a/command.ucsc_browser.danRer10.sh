cwd=`pwd`
datadir=`echo $cwd | sed 's/data_cmd/data/'`

# UCSC genome browser files:
# .2bit

module load ucsc bedtools

cd $datadir

asm=carAur03
asm_dir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/$asm
ucsc_asm_dir=~/data/datashare/fishCarAur1/$asm

asm1=carAur01
asm_dir1=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01

kent=/data/genome/jksrc_v352/kent/src/hg/lib/
chain=$asm2_dir/carAur01_to_$asm2.liftOver.chain
mkdir $asm_dir/big
mkdir $asm_dir/kent_lib
for f in bigChain bigLink bigPsl bigMaf
do
    cp $kent/$f.as $asm_dir/kent_lib/
done

cp ~/data/danRer10/danRer10.chromSizes $asm_dir/danRer10.chromSizes
cp ~/data/common_carp/ng/carp.ng.chromSizes $asm_dir/carp.ng.chromSizes
cut -f1 ~/data/danRer10/danRer10.chromSizes | cut -f1 > $asm_dir/danRer10.chrom.lst
cut -f1 $asm_dir/$asm.chromSizes | cut -f1 > $asm_dir/$asm.chrom.lst
cut -f1 ~/data/common_carp/ng/carp.ng.chromSizes | cut -f1 > $asm_dir/carp.ng.chrom.lst

danRer10_2bit=~/data/danRer10/danRer10.2bit
danRer10_chromsizes=~/data/danRer10/danRer10.chromSizes
carp_ng_2bit=~/data/common_carp/ng/carp.ng.renameLG.2bit
carp_ng_chromsizes=~/data/common_carp/ng/carp.ng.renameLG.fasta.chromSizes

# 
# target=ZF  query=GF  net.chain
in_chain=~/data/goldfish/11549472/sergey_canu70x/arrow/chain_net1/carAur03.vs.danRer10/all.query.syn.net.chain
hgLoadChain -noBin -test $asm bigChain $in_chain
sed 's/.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $asm_dir/danRer10.$asm.net.bigchain
bedToBigBed -type=bed6+6 -as=$asm_dir/kent_lib/bigChain.as -tab $asm_dir/danRer10.$asm.net.bigchain ~/data/danRer10/danRer10.chromSizes $asm_dir/bbi/danRer10.$asm.net.bigchain.bb
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > danRer10.carAur03.net.bigLink
bedToBigBed -type=bed4+1 -as=$asm_dir/kent_lib/bigLink.as -tab danRer10.carAur03.net.bigLink ~/data/danRer10/danRer10.chromSizes bbi/danRer10.carAur03.net.bigLink.bb
# chain to psl
#chainToPsl $in_chain $asm_dir/danRer10.chromSizes $asm_dir/$asm.chromSizes ~/data/danRer10/danRer10.2bit $asm_dir/$asm.2bit > $asm_dir/big/danRer10.carAur03.net.psl
#chainFilter -notQ=tig00037772_arrow $in_chain | chainToAxt -maxGap=100 stdin ~/data/danRer10/danRer10.2bit $asm_dir/$asm.2bit stdout | axtSort stdin stdout | axtToMaf -tPrefix=danRer10. -qPrefix=carAur03. stdin $danRer10_chromsizes $asm_dir/$asm.chromSizes $asm_dir/big/danRer10.carAur03.net.maf;
#mafToBigMaf danRer10 $asm_dir/big/danRer10.carAur03.net.maf stdout | sort -k1,1 -k2,2n > $asm_dir/big/danRer10.carAur03.net.bigMaf
chainFilter -notQ=tig00037772_arrow $in_chain | chainToAxt -maxGap=100 stdin ~/data/danRer10/danRer10.2bit $asm_dir/$asm.2bit stdout | axtToChain stdin $danRer10_chromsizes $asm_dir/$asm.chromSizes stdout | chainToPslBasic stdin stdout | pslToBigPsl -cds=carAur03.CDS.for_bigpsl -fa=$asm_dir/$asm.sm.fa stdin stdout | sort -k1,1 -k2,2n -S 1G > $asm_dir/big/danRer10.carAur03.net.bigPsl 
bedToBigBed -as=$kent/bigPsl.as -type=bed12+13 -tab $asm_dir/big/danRer10.carAur03.net.bigPsl $danRer10_chromsizes $asm_dir/bbi/danRer10.carAur03.net.bigPsl.bb

# target=ZF  query=CC  net.chain
in_chain=~/data/goldfish/11549472/sergey_canu70x/arrow/chain_net1/carp_ng.vs.danRer10/all.query.syn.net.chain
hgLoadChain -noBin -test danRer10 bigChain $in_chain
sed 's/.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $asm_dir/danRer10.carp_ng.net.bigchain
bedToBigBed -type=bed6+6 -as=$asm_dir/kent_lib/bigChain.as -tab $asm_dir/danRer10.carp_ng.net.bigchain ~/data/danRer10/danRer10.chromSizes $asm_dir/bbi/danRer10.carp_ng.net.bigchain.bb
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > danRer10.carp_ng.net.bigLink
bedToBigBed -type=bed4+1 -as=$asm_dir/kent_lib/bigLink.as -tab danRer10.carp_ng.net.bigLink ~/data/danRer10/danRer10.chromSizes bbi/danRer10.carp_ng.net.bigLink.bb
# chain to psl, break at gap > 100bp
chainToAxt -maxGap=100 $in_chain $danRer10_2bit $carp_ng_2bit stdout | \
    axtToChain stdin $danRer10_chromsizes $carp_ng_chromsizes stdout | \
    chainToPslBasic stdin stdout | pslToBigPsl stdin stdout | \
    sort -k1,1 -k2,2n -S 1G > $asm_dir/big/danRer10.carp_ng.net.bigPsl 
bedToBigBed -as=$kent/bigPsl.as -type=bed12+13 -tab $asm_dir/big/danRer10.carp_ng.net.bigPsl $danRer10_chromsizes $asm_dir/bbi/danRer10.carp_ng.net.bigPsl.bb

# target=ZF  query=GF  chain
in_chain=~/data/goldfish/11549472/sergey_canu70x/arrow/chain_net1/carAur03.vs.danRer10/all.swap.chain
hgLoadChain -noBin -test $asm bigChain $in_chain
sed 's/.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $asm_dir/danRer10.$asm.bigchain
bedToBigBed -type=bed6+6 -as=$asm_dir/kent_lib/bigChain.as -tab $asm_dir/danRer10.$asm.bigchain ~/data/danRer10/danRer10.chromSizes $asm_dir/bbi/danRer10.$asm.bigchain.bb
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > danRer10.carAur03.bigLink
bedToBigBed -type=bed4+1 -as=kent_lib/bigLink.as -tab danRer10.carAur03.bigLink ~/data/danRer10/danRer10.chromSizes bbi/danRer10.carAur03.bigLink.bb

in_chain=~/data/goldfish/11549472/sergey_canu70x/arrow/chain_net1/carAur03.vs.danRer10/all.AntiRepeat.swap.chain
hgLoadChain -noBin -test $asm bigChain $in_chain
sed 's/.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $asm_dir/big/danRer10.$asm.AntiRepeat.bigchain
bedToBigBed -type=bed6+6 -as=$kent/bigChain.as -tab $asm_dir/big/danRer10.$asm.AntiRepeat.bigchain $danRer10_chromsizes $asm_dir/bbi/danRer10.$asm.AntiRepeat.bigchain.bb
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > $asm_dir/big/danRer10.carAur03.AntiRepeat.bigLink
bedToBigBed -type=bed4+1 -as=kent_lib/bigLink.as -tab $asm_dir/big/danRer10.carAur03.AntiRepeat.bigLink ~/data/danRer10/danRer10.chromSizes $asm_dir/bbi/danRer10.carAur03.AntiRepeat.bigLink.bb

# chain to psl, break at gap > 100bp
chainToAxt -maxGap=100 $in_chain $danRer10_2bit $carp_ng_2bit stdout | \
    axtToChain stdin $danRer10_chromsizes $carp_ng_chromsizes stdout | \
    chainToPslBasic stdin stdout | pslToBigPsl stdin stdout | \
    sort -k1,1 -k2,2n -S 1G > $asm_dir/big/danRer10.carp_ng.net.bigPsl 

# target=ZF  query=CC  chain
in_chain=~/data/goldfish/11549472/sergey_canu70x/arrow/chain_net1/carp_ng.vs.danRer10/all.swap.chain
hgLoadChain -noBin -test $asm bigChain $in_chain
sed 's/.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > $asm_dir/danRer10.carp_ng.bigchain
bedToBigBed -type=bed6+6 -as=$asm_dir/kent_lib/bigChain.as -tab $asm_dir/danRer10.carp_ng.bigchain ~/data/danRer10/danRer10.chromSizes $asm_dir/bbi/danRer10.carp_ng.bigchain.bb
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > danRer10.carp_ng.bigLink
bedToBigBed -type=bed4+1 -as=kent_lib/bigLink.as -tab danRer10.carp_ng.bigLink ~/data/danRer10/danRer10.chromSizes bbi/danRer10.carp_ng.bigLink.bb


# copy to ucsc dir
ucsc_danRer10_dir=~/data/datashare/fishCarAur1/danRer10
files=""
files=$files bbi/danRer10.carAur03.net.bigchain.bb  bbi/danRer10.carAur03.bigchain.bb
files=$files bbi/danRer10.carAur03.net.bigLink.bb  bbi/danRer10.carAur03.bigLink.bb
files=$files bbi/danRer10.carAur03.AntiRepeat.bigchain.bb  bbi/danRer10.carAur03.AntiRepeat.bigLink.bb  
files=$files bbi/danRer10.carAur03.net.bigPsl.bb bbi/danRer10.carp_ng.bigPsl.bb
files=$files bbi/danRer10.carp_ng.net.bigchain.bb  bbi/danRer10.carp_ng.bigchain.bb
files=$files bbi/danRer10.carp_ng.net.bigLink.bb  bbi/danRer10.carp_ng.bigLink.bb
for f in 
do
    cp $asm_dir/$f $ucsc_danRer10_dir/$f
done
