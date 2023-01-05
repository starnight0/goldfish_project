mkdir big/WGD
cp ../chain_net1/carAur03.vs.carAur03/all.target.syn.net.no_het.chain  big/WGD/GF.GF.net.chain
# WGD1
#cp ../chain_net1/carAur03.vs.carp_ncbi/all.target.syn.net.comb.no_het.chain  big/WGD/GF.CC.net.chain
#cp ../chain_net1/carAur03.vs.grass_carp_female/all.target.syn.net.comb.no_het.chain  big/WGD/GF.GC.net.chain
#cp ../chain_net1/carAur03.vs.danRer10/all.target.syn.net.comb.no_het.chain  big/WGD/GF.ZF.net.chain
# WGD
cp ../chain_net1/carAur03.vs.carp_ncbi/all.target.syn.net.no_het.chain  big/WGD/GF.CC.net.chain
cp ../chain_net1/carAur03.vs.grass_carp_female/all.target.syn.net.no_het.chain  big/WGD/GF.GC.net.chain
cp ../chain_net1/carAur03.vs.danRer10/all.target.syn.net.no_het.chain  big/WGD/GF.ZF.net.chain
cp ../chain_net1/carp_ncbi.vs.danRer10/all.target.syn.net.chain  big/WGD/CC.ZF.net.chain

chainFilter -qMinSize=50000 -tMinSize=50000 big/WGD/GF.GF.net.chain | perl -ne '
if (m/^#/) { next; } chomp; 
if (m/^chain/) {
    my @t=split / /;
    if ($t[2] eq $t[7]) { $filt=1;}
    else {$filt=0; print $_, "\n";}
} elsif (!$filt) {
    print $_, "\n";
}' > big/WGD/GF.GF.net.long.chain
~/my_program3/src/utility/czl_chain_break_at_gap.pl -tg 20000 -qg 20000 -i big/WGD/GF.GF.net.long.chain -o - | chainFilter -qMinSize=50000 -tMinSize=50000 stdin > big/WGD/GF.GF.net.long.G20000.chain
~/my_program3/src/utility/czl_chain_remove_synmetry.pl -ovlf 0.5 -i big/WGD/GF.GF.net.long.G20000.chain -o big/WGD/GF.GF.net.long.G20000.1.chain
~/my_program3/src/utility/czl_chain_remove_TQovl.pl -m TorQ -ovlf 0.5 -i big/WGD/GF.GF.net.long.G20000.1.chain -o big/WGD/GF.GF.net.long.G20000.2.chain
chainSwap big/WGD/GF.GF.net.long.G20000.2.chain big/WGD/GF.GF.net.long.G20000.2.swap.chain
# output big/WGD/GF.GF.net.long.G20000.3.chain with symmetric alignment
chainMergeSort big/WGD/GF.GF.net.long.G20000.2.chain big/WGD/GF.GF.net.long.G20000.2.swap.chain > big/WGD/GF.GF.net.long.G20000.3.chain

chainFilter -qMinSize=50000 -tMinSize=50000 big/WGD/GF.ZF.net.chain > big/WGD/GF.ZF.net.long.chain 
~/my_program3/src/utility/czl_chain_remove_TQovl.pl -m T -ovlf 0.2 -i big/WGD/GF.ZF.net.long.chain -o big/WGD/GF.ZF.net.long.non_Tovl.chain

chainFilter -qMinSize=50000 -tMinSize=50000 big/WGD/CC.ZF.net.chain > big/WGD/CC.ZF.net.long.chain 
~/my_program3/src/utility/czl_chain_remove_TQovl.pl -m T -ovlf 0.2 -i big/WGD/CC.ZF.net.long.chain -o big/WGD/CC.ZF.net.long.non_Tovl.chain

# output bed format
cat big/WGD/GF.GF.net.long.G20000.2.chain | awk -F' ' -v OFS=$'\t' '/^chain/ {print $3,$6,$7,$13,".","+"; if ($10=="+") {print $8,$11,$12,$13,".","+"} else {print $8,$9-$12,$9-$11,$13,".","+"} }' | sort -k1,1 -k2,2n > big/WGD/GF.GF.net.long.G20000.2.chain.bed
# keep only target < query
#cat big/WGD/GF.GF.net.long.G20000.chain | perl -ne '
#if (m/^#/) { next; } chomp; 
#if (m/^chain/) {
#    my @t=split / /;
#    if ($t[2] ge $t[7]) { $filt=1; print $_, "\n";}
#    else {$filt=0;}
#} elsif (!$filt) {
#    print $_, "\n";
#}' > big/WGD/GF.GF.net.long.G20000.half.chain 
#~/my_program3/src/utility/czl_chain_remove_TQovl.pl  -i big/WGD/GF.GF.net.long.G20000.half.chain -o big/WGD/GF.GF.net.long.G20000.half.1.chain -ovlf 0.1 -m TorQ
#chainSwap big/WGD/GF.GF.net.long.G20000.half.1.chain big/WGD/GF.GF.net.long.G20000.half.1.swap.chain

#############################
# combine all exon/CNE annotation into a bed file
#############################
# {{{
for sp in ZF CC GF
do
	if [ "$sp" == "ZF" ]
	then
		cat ~/data/zebrafish/ensembl85/Danio_rerio.GRCz10.85.gtf.exon.gte.bed | sort -k1,1 -k2,2n > big/WGD/$sp.exon.bed
		cat big/phastcons/run2.net_roast_multic/CC_ZF.ZF.notRE.bed big/phastcons/run2.net_roast_multic/CC_ZF.ZF.notRE.bed
	fi
	if [ "$sp" == "CC" ]
	then
		cat ~/data/common_carp/NCBI_Cyprinus_carpio/ref_common_carp_genome_top_level.rename_chr.gff3.exon.gte.bed | sort -k1,1 -k2,2n > big/WGD/$sp.exon.bed
	fi
	if [ "$sp" == "GF" ]
	then
		cat big/carAur03.noM.exon.bed  big/carAur03.noM.ncrna.f.not_ovl.exon.bed | sort -k1,1 -k2,2n | sed 's/\t\(CA.*\)_R/\t\1:\1_R/' > big/WGD/$sp.exon.bed
	fi
	cat big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.$sp.notR.bed | sort -k1,1 -k2,2n > big/WGD/$sp.CNE.notR.bed
	bedtools subtract -nonamecheck -A -a big/WGD/$sp.CNE.notR.bed -b big/WGD/$sp.exon.bed > big/WGD/$sp.CNE.bed 

#cat big/carAur03.noM.exon.bed  big/carAur03.noM.ncrna.f.not_ovl.exon.bed big/WGD/CNE.bed | bedtools subtract -nonamecheck -A -a big/phastcons/run2.net_roast_multic/GF_GF.notRE.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.GF_GF.bed
cat $sp.exon.bed $sp.CNE.bed | bedtools subtract -nonamecheck -A -a big/phastcons/run2.net_roast_multic/CC_GF.GF.notRE.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.CC_GF.bed
cat big/carAur03.noM.exon.bed  big/carAur03.noM.ncrna.f.not_ovl.exon.bed big/WGD/CNE.bed | bedtools subtract -nonamecheck -A -a big/phastcons/run2.net_roast_multic/GC_GF.GF.notRE.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.GC_GF.bed
cat big/carAur03.noM.exon.bed  big/carAur03.noM.ncrna.f.not_ovl.exon.bed big/WGD/CNE.bed | bedtools subtract -nonamecheck -A -a big/phastcons/run2.net_roast_multic/ZF_GF.GF.notRE.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.ZF_GF.bed

#cat big/WGD/CNE.CC_GF.bed big/WGD/CNE.ZF_GF.bed | bedtools subtract -nonamecheck -A -a big/WGD/CNE.GF_GF.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.GF_GF_only.bed
cat big/WGD/CNE.ZF_GF.bed big/WGD/CNE.GC_GF.bed | bedtools subtract -nonamecheck -A -a big/WGD/CNE.CC_GF.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.CC_GF.1.bed
cat big/WGD/CNE.CC_GF.1.bed big/WGD/CNE.ZF_GF.bed | bedtools subtract -nonamecheck -A -a big/WGD/CNE.GC_GF.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.GC_GF.1.bed
cat big/WGD/CNE.CC_GF.1.bed big/WGD/CNE.GC_GF.1.bed | bedtools subtract -nonamecheck -A -a big/WGD/CNE.ZF_GF.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.ZF_GF.1.bed

>big/WGD/exon_CNE.bed.tmp;
cat  big/WGD/exon.bed | awk -F$'\t' -v OFS=$'\t' '{gsub(/:/, "__",$4); $4="exon:"$4":"$1":"$2":"$3; print}' > big/WGD/exon_CNE.bed.tmp;
cat  big/WGD/CNE.bed  | awk -F$'\t' -v OFS=$'\t' '{$4="CNE1:"$4":"$1":"$2":"$3; print}' >> big/WGD/exon_CNE.bed.tmp
# cat  big/WGD/CNE.GF_GF_only.bed  | awk -F$'\t' -v OFS=$'\t' '{$4="CNE2:"$4":"$1":"$2":"$3; print}' >> big/WGD/exon_CNE.bed.tmp
cat  big/WGD/CNE.CC_GF.1.bed  | awk -F$'\t' -v OFS=$'\t' '{$4="CNE3:"$4":"$1":"$2":"$3; print}' >> big/WGD/exon_CNE.bed.tmp
cat  big/WGD/CNE.GC_GF.1.bed  | awk -F$'\t' -v OFS=$'\t' '{$4="CNE4:"$4":"$1":"$2":"$3; print}' >> big/WGD/exon_CNE.bed.tmp
cat  big/WGD/CNE.ZF_GF.1.bed  | awk -F$'\t' -v OFS=$'\t' '{$4="CNE5:"$4":"$1":"$2":"$3; print}' >> big/WGD/exon_CNE.bed.tmp
sort -k1,1 -k2,2n big/WGD/exon_CNE.bed.tmp >  big/WGD/exon_CNE.bed
rm big/WGD/exon_CNE.bed.tmp
join -a 1 -t$'\t' -j1 big/WGD/exon_CNE.bed ../carAur01/carAur01.masked_ctg_id  | awk -F$'\t' 'NF==6 && $1!="chrM"' > big/WGD/exon_CNE.bed.tmp
mv big/WGD/exon_CNE.bed.tmp big/WGD/exon_CNE.bed 

bedtools intersect -nonamecheck -a big/WGD/exon_CNE.bed -b ../JASPAR_run1/fimo_out/carAur03.fimo.1.bed -wao | awk -F$'\t' -v OFS=$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6,"."} $7!="." { if (chr==$1 && b==$2 && id==$4) {id2=id2","$10} else {if (id2) {print chr,b,e,id,sc,s,id2} chr=$1; b=$2; e=$3; id=$4; sc=$5; s=$6; id2=$10; } } END {if (id2) {print chr,b,e,id,sc,s,id2}}' > big/WGD/exon_CNE.add_motif.bed

# }}}
#################################

cat big/WGD/exon_CNE.bed | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.bed.tmp

crossmap bed big/WGD/GF.GF.net.chain big/WGD/exon_CNE.bed  big/WGD/exon_CNE.map.bed
cat big/WGD/exon_CNE.map.bed | sort -k1,1 -k2,2n -k4,4 | awk -F$'\t' -v OFS=$'\t' '{if (id!=$4 || chr!=$1) {if (a) {print a;}  a=$0; id=$4;chr=$1; } } END{ if (a) {print a} }' | sort -k4,4 -k1,1 -k2,2n | awk -F$'\t' -v OFS=$'\t' '{if (id!=$4 || chr!=$1 || $2-e>1000 ) {if (a) {print a;}  a=$0; id=$4;chr=$1;e=$3; } else {e=$3;} } END{ if (a) {print a} }' > big/WGD/exon_CNE.map.1.bed
bedtools intersect -nonamecheck -a big/WGD/exon_CNE.map.1.bed -b big/WGD/exon_CNE.bed  -wao | awk -F$'\t' -v OFS=$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6,"."} $7!="." { if (chr==$1 && b==$2 && id==$4) {id2=id2","$10} else {if (id2) {print chr,b,e,id,sc,s,id2} chr=$1; b=$2; e=$3; id=$4; sc=$5; s=$6; id2=$10; } } END {if (id2) {print chr,b,e,id,sc,s,id2}}' > big/WGD/exon_CNE.map.2.bed
cat big/WGD/exon_CNE.map.2.bed | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.map.bed.tmp
join -a 1 -e"." -t$'\t' -j1  -o 1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.8,2.6,2.7  big/WGD/exon_CNE.bed.tmp   big/WGD/exon_CNE.map.bed.tmp  > big/WGD/exon_CNE.AB.bed

for sp in CC ZF GC
do
crossmap bed big/WGD/GF.$sp.net.chain big/WGD/exon_CNE.bed  big/WGD/exon_CNE.map.$sp.bed
cat big/WGD/exon_CNE.map.$sp.bed | sort -k1,1 -k2,2n -k4,4 | awk -F$'\t' -v OFS=$'\t' '{if (id!=$4 || chr!=$1) {if (a) {print a;}  a=$0; id=$4;chr=$1; } } END{ if (a) {print a} }' | sort -k4,4 -k1,1 -k2,2n | awk -F$'\t' -v OFS=$'\t' '{if (id!=$4 || chr!=$1 || $2-e>1000 ) {if (a) {print a;}  a=$0; id=$4;chr=$1;e=$3; } else {e=$3;} } END{ if (a) {print a} }' > big/WGD/exon_CNE.map.$sp.1.bed
cat big/WGD/exon_CNE.map.$sp.1.bed | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.map.$sp.bed.tmp
join -a 1 -e"." -t$'\t' -j1  -o 1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.5,2.6,2.7  big/WGD/exon_CNE.bed.tmp   big/WGD/exon_CNE.map.$sp.bed.tmp  > big/WGD/exon_CNE.map.$sp.bed12
done


#########################
# split chain to each chain per file
#########################
mkdir big/WGD/GF.GF.net.long.G20000.chain.split
mkdir big/WGD/exon.split
mkdir big/WGD/CNE.split
mkdir big/WGD/exon_CNE.split
cat big/WGD/GF.GF.net.long.G20000.3.chain | awk -v i=1 '/^chain/ || /^[0-9]/ {print >> "big/WGD/GF.GF.net.long.G20000.chain.split/"i".chain"} /^\s*$/ {print $0 >> "big/WGD/GF.GF.net.long.G20000.chain.split/"i".chain"; i++; }'
cat big/WGD/GF.GF.net.long.G20000.3.chain | grep 'chain' | awk -v OFS=$'\t' '{print $3,$6,$7}' > big/WGD/GF.GF.net.long.G20000.bed3
cat big/WGD/GF.GF.net.long.G20000.3.chain | grep 'chain' | awk -v OFS=$'\t' '$10=="+" {print $8,$11,$12} $10=="-" {print $8,$9-$12,$9-$11}' > big/WGD/GF.GF.net.long.G20000.Q.bed3
###########################
# liftover per chain
###########################
n=`cat big/WGD/GF.GF.net.long.G20000.bed3 | wc -l`
IFS=$'\n'
for i in `seq 1 $n`
do
    if ! [ -f big/WGD/GF.GF.net.long.G20000.chain.split/$i.chain ]; then continue; fi;

    echo $i;
    chr=`cut -f 1 big/WGD/GF.GF.net.long.G20000.bed3 | head -n $i | tail -n 1`;
    b=`cut -f 2 big/WGD/GF.GF.net.long.G20000.bed3 | head -n $i | tail -n 1`;
    e=`cut -f 3 big/WGD/GF.GF.net.long.G20000.bed3 | head -n $i | tail -n 1`;
    qchr=`cut -f 1 big/WGD/GF.GF.net.long.G20000.Q.bed3 | head -n $i | tail -n 1`;
    qb=`cut -f 2 big/WGD/GF.GF.net.long.G20000.Q.bed3 | head -n $i | tail -n 1`;
    qe=`cut -f 3 big/WGD/GF.GF.net.long.G20000.Q.bed3 | head -n $i | tail -n 1`;
#    cat big/WGD/exon.f.bed | awk '$1=="'$chr'" && $2>='$b' && $3<='$e'' | sort -k1,1 -k2,2n > big/WGD/exon.split/$i.A.bed
#    cat big/WGD/CNE.f.bed  | awk '$1=="'$chr'" && $2>='$b' && $3<='$e'' | sort -k1,1 -k2,2n > big/WGD/CNE.split/$i.A.bed
#    cat  big/WGD/exon.split/$i.A.bed | awk -F$'\t' -v OFS=$'\t' '{$4="exon:"$4; print}' > big/WGD/exon_CNE.split/$i.A.bed.tmp
#    cat  big/WGD/CNE.split/$i.A.bed  | awk -F$'\t' -v OFS=$'\t' '{$4="CNE:"$4; print}' >> big/WGD/exon_CNE.split/$i.A.bed.tmp
    cat big/WGD/exon_CNE.bed | awk -F$'\t' '$1=="'$chr'" && $2>='$b' && $3<='$e'' > big/WGD/exon_CNE.split/$i.A.bed
    cat big/WGD/exon_CNE.bed | awk -F$'\t' '$1=="'$qchr'" && $2>='$qb' && $3<='$qe'' > big/WGD/exon_CNE.split/$i.B.bed
    crossmap bed  big/WGD/GF.GF.net.long.G20000.chain.split/$i.chain  big/WGD/exon_CNE.split/$i.A.bed  big/WGD/exon_CNE.split/$i.A.map.bed
    cat big/WGD/exon_CNE.split/$i.A.map.bed | sort -k1,1 -k2,2n -k4,4 | awk -F$'\t' -v OFS=$'\t' '{ if (chr!=$1 || id!=$4) { if (a) {print a;} chr=$1; id=$4; a=$0;}  } END{ if (a) {print a;} }' | sort -k4,4 -k1,1 -k2,2n | awk -F$'\t' -v OFS=$'\t' '{if (id!=$4 || chr!=$1 || $2-e>1000 ) {if (a) {print a;}  a=$0; id=$4;chr=$1;e=$3; } else {e=$3;} } END{ if (a) {print a} }' > big/WGD/exon_CNE.split/$i.A.map.1.bed;
    bedtools intersect -nonamecheck -a big/WGD/exon_CNE.split/$i.A.map.1.bed -b big/WGD/exon_CNE.split/$i.B.bed -wao | awk -F$'\t' -v OFS=$'\t' '$7=="." {print $1,$2,$3,$4,$5,$6,"."} $7!="." { if (chr==$1 && b==$2 && id==$4) {id2=id2","$10} else {if (id2) {print chr,b,e,id,sc,s,id2} chr=$1; b=$2; e=$3; id=$4; sc=$5; s=$6; id2=$10; } } END {if (id2) {print chr,b,e,id,sc,s,id2}}' > big/WGD/exon_CNE.split/$i.A.map.2.bed
    cat big/WGD/exon_CNE.split/$i.A.map.2.bed | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 -k3,3n > big/WGD/exon_CNE.split/$i.A.map.bed.tmp
    cat big/WGD/exon_CNE.split/$i.A.bed     | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 -k3,3n > big/WGD/exon_CNE.split/$i.A.bed.tmp
    join -a 1 -e"." -t$'\t' -j1  -o 1.2,1.3,1.4,1.5,1.6,1.7,2.2,2.3,2.4,2.8,2.6,2.7  big/WGD/exon_CNE.split/$i.A.bed.tmp   big/WGD/exon_CNE.split/$i.A.map.bed.tmp  > big/WGD/exon_CNE.split/$i.AB.bed
    rm big/WGD/exon_CNE.split/$i.A.map.bed.tmp
    rm big/WGD/exon_CNE.split/$i.A.bed.tmp

done

>big/WGD/exon_CNE.by_chain.AB.txt
>big/WGD/exon_CNE.by_chain.AB.counts.txt
for i in `seq 1 $n`
do
    cat big/WGD/exon_CNE.split/$i.AB.bed | awk -F$'\t' '{print '$i'"\t"$0}' >> big/WGD/exon_CNE.by_chain.AB.txt;
    m1=`cat big/WGD/exon_CNE.split/$i.AB.bed | cut -f 1-6 | sort | uniq | wc -l`;
    m2=`cat big/WGD/exon_CNE.split/$i.AB.bed | awk -F$'\t' '$7=="."' | cut -f 1-6 | sort | uniq | wc -l`;
    echo $m1$'\t'$m2 >> big/WGD/exon_CNE.by_chain.AB.counts.txt;
done


cat big/WGD/exon_CNE.by_chain.AB.txt | awk -F$'\t' -v OFS=$'\t' '{print $5"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.by_chain.AB.txt.tmp
cat big/WGD/exon_CNE.by_chain.AB.txt | awk -F$'\t' -v OFS=$'\t' '$8=="." {print $5"\t"$1,$2,$3,$4,$5,$6,$7}' | sort -k1,1 > big/WGD/exon_CNE.by_chain.AB.unmap.txt.tmp
cat big/WGD/exon_CNE.AB.bed | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.AB.bed.tmp
join -t$'\t' -j1  -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.8,2.9,2.10,2.11,2.12,2.13  big/WGD/exon_CNE.by_chain.AB.unmap.txt.tmp big/WGD/exon_CNE.AB.bed.tmp | sort -k1,1 > big/WGD/a1.tmp
join -a1 -t$'\t' -j1 -e "." -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.9,2.10,2.11,2.12,2.13,2.14   big/WGD/exon_CNE.by_chain.AB.txt.tmp   big/WGD/a1.tmp | sort -k1,1n -k2,2 -k3,3n  > big/WGD/exon_CNE.by_chain.AB.tab19.txt

cat big/WGD/exon_CNE.by_chain.AB.tab19.txt | awk -F$'\t' '{print $5"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.by_chain.AB.tab19.txt.tmp
cat big/WGD/exon_CNE.map.CC.bed12 | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.map.CC.bed12.tmp
cat big/WGD/exon_CNE.map.GC.bed12 | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.map.GC.bed12.tmp
cat big/WGD/exon_CNE.map.ZF.bed12 | awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 > big/WGD/exon_CNE.map.ZF.bed12.tmp
join -t$'\t' -j1 big/WGD/exon_CNE.by_chain.AB.tab19.txt.tmp  big/WGD/exon_CNE.map.CC.bed12.tmp | cut -f 1-20,27-32 > big/WGD/exon_CNE.by_chain.AB.tab25.txt.tmp
cut -f 2- big/WGD/exon_CNE.by_chain.AB.tab25.txt.tmp | sort -k1,1n -k2,2 -k3,3n > big/WGD/exon_CNE.by_chain.AB.tab25.txt
join -t$'\t' -j1 big/WGD/exon_CNE.by_chain.AB.tab25.txt.tmp  big/WGD/exon_CNE.map.GC.bed12.tmp | cut -f 1-26,33-38 > big/WGD/exon_CNE.by_chain.AB.tab31.txt.tmp
cut -f 2- big/WGD/exon_CNE.by_chain.AB.tab31.txt.tmp | sort -k1,1n -k2,2 -k3,3n > big/WGD/exon_CNE.by_chain.AB.tab31.txt
join -t$'\t' -j1 big/WGD/exon_CNE.by_chain.AB.tab31.txt.tmp  big/WGD/exon_CNE.map.ZF.bed12.tmp | cut -f 1-32,39-44 > big/WGD/exon_CNE.by_chain.AB.tab37.txt.tmp
cut -f 2- big/WGD/exon_CNE.by_chain.AB.tab37.txt.tmp | sort -k1,1n -k2,2 -k3,3n > big/WGD/exon_CNE.by_chain.AB.tab37.txt
cat big/WGD/exon_CNE.by_chain.AB.tab37.txt | awk '$8!="." || $14!="." || $20!="." || $26!="." || $32!="." || $5~/exon/' > big/WGD/exon_CNE.by_chain.AB.tab37.1.txt
cat big/WGD/exon_CNE.by_chain.AB.tab37.1.txt | sort -k2,2 -k3,3n -k5,5 -k1,1n > big/WGD/exon_CNE.by_chain.AB.tab37.2.txt 


for len in 2 5 10 20
do
cat big/WGD/exon_CNE.by_chain.AB.tab37.2.txt | perl -e 'my %v; my %count; my %pass;
while(<>) {chomp; my @t=split "\t"; my $aid=$t[0]; push @{$v{$aid}}, \@t;}
my $m = 5;
foreach my $aid (keys(%v)) {
	my $v=$v{$aid};
for (my $i=0; $i<@$v; $i++) {
    if ($v->[$i][7] eq "." && $v->[$i][13] eq ".") {
        my ($aid, $chr, $b, $e, $id) = @{$v->[$i]}[0..4];
        if ($i>0 && $v->[$i-1][7] eq "." && $v->[$i-1][1] eq $chr && $v->[$i-1][2]==$b && $v->[$i-1][4] eq $id) { next; }
		if (exists $pass{$id}) {next; }
		$pass{$id}++;
        my $i0 = $i;   my $gid0;  my $d0;
        my $i1 = $i;   my $gid1;  my $d1;
        if ($id=~m/^exon/i) {
            if ($id=~m/exon:([^:]+)_R/) { $gid0=$1; $d0=0; }
        } else {
            while ($i0>0) {
                $i0--;
                if ($v->[$i0][4]=~m/exon:([^:]+)_R/) {
                    $gid0=$1; $d0=$b-$v->[$i0][3]; last;
                }
            }
            while ($i1<@v) {
                if ($v->[$i1][4]=~m/exon:([^:]+)_R/) {
                    $gid1=$1; $d1=$v->[$i1][2]-$3; last;
                }
                $i1++;
            }
            if (defined $gid0 && defined $gid1 && $gid0 eq $gid1) { undef $gid1; undef $d1; $d0=0;}
            else {
                if (defined $gid0 && $d0>'$len'000) { undef $gid0; undef $d0; }
                if (defined $gid1 && $d1>'$len'000) { undef $gid1; undef $d1; }
            }
        }
        my @k;
        my $k0;
        if ($id=~/exon/) { $k0 = 0; }
        else { $k0 = 4*$m; }
        push @k, $k0;
        if ($v->[$i][19] eq "." && $v->[$i][25] eq ".") {
            if ($v->[$i][31] eq ".") { push @k,$k0+1; }
            else { push @k, $k0+3; }
        } else {
            if ($v->[$i][31] eq ".") { push @k,$k0+2; }
            else { push @k, $k0+4; }
        }

        foreach my $k (@k) {
            if (defined $gid0) {
                my $w0;
                if ($d0==0) { $w0=2; }
                elsif ($d0<1000) { $w0=2-$d0/1000; }
				else { $w0=1000/$d0; }
                $count{$gid0}[$k]++;
                $count{$gid0}[$m+$k]+=$e-$b;
                $count{$gid0}[$k+$m*2]+=$w0;
                push @{$count{$gid0}[$k+$m*3]}, $id;
            }
            if (defined $gid1) {
                my $w1;
                if ($d1==0) { $w1=2; }
                elsif ($d1<1000) { $w1=2-$d1/1000; }
				else { $w1=1000/$d1; }
                $count{$gid1}[$k]++;
                $count{$gid1}[$k+$m]+=$e-$b;
                $count{$gid1}[$k+$m*2]+=$w1;
                push @{$count{$gid1}[$k+$m*3]}, $id;
            }
        }
    }
}
}
foreach my $gid (sort keys(%count)) {
    print $gid;
    for (my $k0=0; $k0<2; $k0++) {
        for (my $k=0; $k<$m*3; $k++) {
            if (!defined $count{$gid}[$k0*$m*4+$k]) { $count{$gid}[$k0*$m*4+$k]=0; }
        }
        print "\t", join("\t", @{$count{$gid}}[($k0*$m*4)..($k0*$m*4+$m*3-1)]);
        for (my $k=$m*3; $k<$m*4; $k++) {
            print "\t", join(",", @{$count{$gid}[$k0*$m*4+$k]});
        }
    }
    print "\n";
}
' > big/WGD/exon_CNE.by_chain.AB.tab37.gene_count.${len}K.txt 
done


cat big/WGD/exon_CNE.by_chain.AB.tab19.txt | sort -k2,2 -k3,3n | perl -e 'my @v; my %count;
while(<>) {chomp; my @t=split "\t"; push @v, \@t;}
for (my $i=0; $i<@v; $i++) {
    if ($v[$i][7] eq ".") {
        my ($aid, $chr, $b, $e, $id) = @{$v[$i]}[0..4];
        if ($i>0 && $v[$i-1][7] eq "." && $v[$i-1][0]==$aid && $v[$i-1][1] eq $chr && $v[$i-1][4] eq $id) { next; }
        my $i0 = $i;   my $gid0;
        while ($i0>0) {
            $i0--;
            if ($v[$i0][0]!=$aid || $b-$v[$i0][3]>20000) { last; }
            if ($v[$i0][4]=~m/exon:([^:]+)_R/) { $gid0=$1; last; }
        }
        my $i1 = $i;   my $gid1;
        while ($i1<@v) {
            if ($v[$i1][0]!=$aid || $v[$i1][2]-$e>20000) { last; }
            if ($v[$i1][4]=~m/exon:([^:]+)_R/) {$gid1=$1; last;}
            $i1++;
        }
        my $is_lost=0;
        if ($v[$i][13] eq ".") { $is_lost=1; }
        if (defined $gid0) { print join("\t", (@{$v[$i]}[1..6], $gid0, $is_lost)), "\n"; }
        if (defined $gid1) { print join("\t", (@{$v[$i]}[1..6], $gid1, $is_lost)), "\n"; }
    }
}
' > big/WGD/exon_CNE.lost.txt

####################################

cat big/WGD/exon_CNE.by_chain.AB.txt | cut -f 5,11 | awk -F$'\t' '$1~/exon/ && $2!="." {print $1}' | awk -F$'\t' '$1~/CA[0-9]/ {gsub(/^exon:/,"",$1); gsub(/:.*$/,"",$1); print $1}' | sed 's/_R[0-9]\+//g' | sort | uniq | wc
cat big/WGD/exon_CNE.by_chain.AB.txt | cut -f 5,11 | awk -F$'\t' '$1~/exon/ && $2~/exon/' | perl -ne 'chomp; my @t=split /\t/; if ($t[1]=~m/,/) { foreach my $a (split /,/,$t[1]) { if ($a=~m/^CNE/i) {next; } print "$t[0]\t$a\n"; last; } } else {print $_,"\n";}' > a1
awk -F$'\t' -v OFS=$'\t' '$1<$2 {print} $1>$2 {print $2,$1}' a1 | sort -k1,1 -k2,2 | uniq > a2;
cat a2 | awk -F$'\t' '$1~/CA[0-9]/ && $2~/CA[0-9]/ {gsub(/^exon:/,"",$1); gsub(/:.*$/,"",$1); gsub(/^exon:/,"",$2); gsub(/:.*$/,"",$2); print $1"\t"$2}' | sed 's/_R[0-9]\+//g' | sort -k1,1 -k2,2 | uniq -c | awk -v OFS=$'\t' '{print $2,$3,$1}' > gene_pairs0.txt
cat gene_pairs0.txt | sort -k1,1 -k3,3nr | awk -v OFS=$'\t' '{ if (id1!=$1) {print; id1=$1;} }' | sort -k2,2 -k3,3nr | awk -v OFS=$'\t' '{ if (id2!=$2) {print $1,$2; id2=$2; } }' > big/WGD/gene_pairs.txt
cat big/carAur03.noM.gene.bgp | cut -f 1,18,19 | awk -F$'\t' '{print $2"\t"$0}' | sort -k1,1 > a3
cat big/WGD/gene_pairs.txt | sort -k1,1 > b1
join -t$'\t' -j1 -o 2.2,1.2,1.3,1.4 a3 b1 | sort -k1,1 > b2
join -t$'\t' -j1 -o 2.2,2.3,2.4,1.2,1.3,1.4 a3 b2 | sort -k1,1 -k2,2 > big/WGD/gene_pairs.1.txt

####################################
cat  big/WGD/exon_CNE.by_chain.AB.tab37.2.txt | perl -e '
my @count = (0) x 12;
my @prev;
my @heads = qw(TotalExon TotalExonLost TotalExonLostNew TotalExonLostMapCC TotalExonLostMapZF TotalExonLostMapCCZF    TotalCNE TotalCNELost TotalCNELostNew TotalCNELostMapCC TotalCNELostMapZF TotalCNELostMapCCZF);
while(<>) {
chomp; my @t=split "\t";
my ($aid, $chr, $b, $e, $id) = @t[0..4];
if ($#prev>=0 && $prev[1] eq $chr && $prev[2]==$b && $prev[4] eq $id) { next; }
if ($t[7] eq "." && $t[13] eq ".") {
    if ($t[4]=~m/exon/) { $k0=0; }
    else { $k0=6;}
    if ($t[19] eq "." && $t[25] eq ".") {
        if ($t[31] eq ".") { $count[$k0+2]++; } else { $count[$k0+4]++; }
    } else {
        if ($t[31] eq ".") { $count[$k0+3]++; } else { $count[$k0+5]++; }
    }
    $count[$k0+1]++;
    $count[$k0]++;
} else { $count[$k0]++; }
@prev = @t;
}
print join("\t", @head  ), "\n";
print join("\t", @count ), "\n"; ' > big/WGD/exon_CNE.by_chain.AB.tab37.count_lost.txt

####################################
# singleton genes from exon_CNE.by_chain.AB.txt
cat big/WGD/exon_CNE.by_chain.AB.txt | perl -e 'my %gc; my %gc1;
while (<>) {
	chomp; @t=split /\t/; my $id1=$t[4]; $id2=$t[10]; if ($id1=~m/^CNE/) {next; }
	if ($id1=~m/^exon:(CA[0-9]+).*$/) { $id1=$1;}
	else { next; }
	if ($id2 eq ".") {
	} else {
		if ($id2=~m/^exon:(CA[0-9]+).*$/) { $id2=$1;}
		$gc1{$id1}++;
	}
	$gc{$id1}++;
}
my $n=0;
foreach my $id (sort keys(%gc)) {
	if (! exists $gc1{$id}) { print "$id\n"; $n++; };
}
print "$n\n";
'

####################################
# number of gene with ZF and retained GF duplicate
####################################
# 60889
cat big/WGD/exon_CNE.by_chain.AB.tab37.2.txt | awk -F$'\t' -v OFS=$'\t' '$32!="." && $2!="." && $8!="." && ($5~/exon/ && $11~/exon/) {split($5,a,":"); split($11,b,":"); print a[2]"\t"b[2];}' | sort | uniq >big/WGD/a2;
wc big/WGD/a2;
cat big/WGD/exon_CNE.by_chain.AB.tab37.2.txt | awk -F$'\t' -v OFS=$'\t' '$32!="." && $2!="." && $8!="." && ($5~/exon/ || $11~/exon/) { if ($5~/exon/) {split($5,a,":"); print a[2];} }' | sort | uniq >big/WGD/a3;
cat big/WGD/exon_CNE.by_chain.AB.tab37.2.txt | awk -F$'\t' -v OFS=$'\t' '$32!="." { if ($2=="." && $11~/exon/) {print $11;} else if ($8=="." && $5~/exon/) {print $5} }' | cut -f 2 -d':' | sort | uniq > big/WGD/a1
cat big/WGD/a1 big/WGD/a3  | sort | uniq -d > big/WGD/a4;
cat big/WGD/a1 big/WGD/a4 | sort | uniq -u > big/WGD/a5; mv big/WGD/a5 big/WGD/a1;

#####################################
# retained WGD gene, singleton gene
#####################################
~/my_program3/src/utility/czl_chain_break_at_gap.pl -tg 20000 -qg 20000 -i big/WGD/GF.GF.net.chain -o big/WGD/GF.GF.net.G20000.chain
~/my_program3/src/utility/czl_chain_remove_synmetry.pl -ovlf 0.5 -i big/WGD/GF.GF.net.G20000.chain -o big/WGD/GF.GF.net.G20000.1.chain
~/my_program3/src/utility/czl_chain_remove_TQovl.pl -m TorQ -ovlf 0.5 -i big/WGD/GF.GF.net.G20000.1.chain -o big/WGD/GF.GF.net.G20000.2.chain
cat big/WGD/GF.GF.net.G20000.2.chain | awk -F' ' -v OFS=$'\t' '/^chain/ {print $3,$6,$7,$13,".","+"; if ($10=="+") {print $8,$11,$12,$13,".","+"} else {print $8,$9-$12,$9-$11,$13,".","+"} }' | sort -k1,1 -k2,2n > big/WGD/GF.GF.net.G20000.2.chain.bed
cat big/carAur03.noM.gene.unmasked.bgp 
bedtools intersect -nonamecheck -a big/carAur03.noM.gene.unmasked.bed -b big/WGD/GF.GF.net.long.G20000.2.chain.bed -wa  -u > big/WGD/GF.unmasked.dup.1.bed
bedtools intersect -nonamecheck -a big/carAur03.noM.gene.bed -b big/WGD/GF.GF.net.long.G20000.2.chain.bed -wa  -u > big/WGD/GF.dup.1.bed
bedtools intersect -nonamecheck -a big/carAur03.noM.gene.unmasked.bed -b big/WGD/GF.GF.net.G20000.2.chain.bed -wa  -u > big/WGD/GF.dup.bed
bedtools subtract -nonamecheck -A -a big/carAur03.noM.gene.unmasked.bed -b big/WGD/GF.GF.net.G20000.2.chain.bed -wa  -u > big/WGD/GF.singleton.bed
