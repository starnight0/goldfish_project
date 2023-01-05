#CWD big/landscape

mkdir big/landscape
mkdir bbi/landscape
mkdir ~/data/datashare/fishCarAur1/carAur01/bbi/landscape
setfacl -m user:webcpu:r-x ~/data/datashare/fishCarAur1/carAur01/bbi/landscape

############# liftOver unmasked contig #############

############# Create various windows ##############
cat carAur03.sm.fa.fai  | awk '{print $1"\t1\t"$2}' | sort -k1,1 -k2,2n > big/landscape/windows.chr.bed
#liftOver ../carAur01/carAur01.unmasked_ctg_id.bed3 $chain carAur03.het_unmasked.bed3 carAur03.het_unmasked.unmap 
CrossMap.py bed $chain ../carAur01/carAur01.unmasked_ctg_id.bed3 big/carAur03.het_unmasked.bed3;
sort -k1,1 -k2,2n big/carAur03.het_unmasked.bed3 > big/carAur03.het_unmasked.bed3.tmp; mv big/carAur03.het_unmasked.bed3.tmp big/carAur03.het_unmasked.bed3;
bedtools subtract -a big/carAur03.het_unmasked.bed3 -b big/carAur03.repeatmasker.merged.bed3 > big/carAur03.het_unmasked.notR.bed3 

bedtools makewindows -g carAur03.sm.fa.fai -w 10000   -s 1000   > big/landscape/windows.10K_1K.bed
bedtools makewindows -g carAur03.sm.fa.fai -w 100000  -s 10000  > big/landscape/windows.100K_10K.bed
bedtools makewindows -g carAur03.sm.fa.fai -w 1000000 -s 10000  > big/landscape/windows.1M_10K.bed
bedtools makewindows -g carAur03.sm.fa.fai -w 1000000 -s 100000 > big/landscape/windows.1M_100K.bed

############# Create track for retained ohnolog regions ###########
#cp ../../../chain_net1/carAur03.vs.carAur03/all.target.syn.net.no_het.chain.20.merged.bed3 ./carAur03.retained_ohnolog.bed
scp $GRYPHON:/data/projects/burgess/zelin/goldfish/11549472/sergey_canu70x/arrow/chain_net1/carAur03.vs.carAur03/all.target.syn.net.no_het.chain.20.merged.bed3 big/landscape/carAur03.retained_ohnolog.bed
bedtools intersect -a big/landscape/carAur03.retained_ohnolog.bed -b big/carAur03.het_unmasked.bed3  > big/landscape/carAur03.retained_ohnolog.no_het.bed 
bedtools subtract -a big/carAur03.het_unmasked.bed3 -b big/landscape/carAur03.retained_ohnolog.no_het.bed > big/landscape/carAur03.lost_ohnolog.no_het.bed;

sort -k1,1 -k2,2n big/landscape/carAur03.retained_ohnolog.bed >big/landscape/carAur03.retained_ohnolog.bed.tmp; mv big/landscape/carAur03.retained_ohnolog.bed.tmp big/landscape/carAur03.retained_ohnolog.bed  
bedToBigBed big/landscape/carAur03.retained_ohnolog.bed carAur03.sm.fa.fai bbi/landscape/carAur03.retained_ohnolog.bb



############ copy files to hub ##########
for f in `ls bbi/landscape/*`; do f1=`echo $f | sed 's/carAur03/carAur01/'`; cp $f ~/data/datashare/fishCarAur1/carAur01/$f1 ; done;
setfacl -m user:webcpu:r-- ~/data/datashare/fishCarAur1/carAur01/bbi/landscape/*
