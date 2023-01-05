## concat ensembl85 fish and GC,GF,CC pep and cdna sequence
>fish17.pep.fa
>fish17.cdna.fa
>fish17.rna_to_gene.map
>fish17.protein_to_gene.map
for i in `seq 1 17`
do
    sp=`head -n $i fish.txt | tail -n 1 | cut -f 4`
    pep=`head -n $i fish.txt | tail -n 1 | cut -f 5`
    cdna=`head -n $i fish.txt | tail -n 1 | cut -f 6`
    gtp=`head -n $i fish.txt | tail -n 1 | cut -f 8`
    zcat $pep  | sed -e '/^>/ s/\s.*$//' -e '/^>/ s/^>/>'$sp'|/' >> fish17.pep.fa
    zcat $cdna | sed -e '/^>/ s/\s.*$//' -e '/^>/ s/^>/>'$sp'|/' >> fish17.cdna.fa
    if [ "$sp" == "GC" ]
    then
        cat $gtp | awk -F$'\t' -v OFS=$'\t' '{gsub(/\..*$/,"",$1); gsub(/\..*$/,"",$2); print "'$sp'|"$2"\t'$sp'|"$1}' >> fish17.rna_to_gene.map
        cat $gtp | awk -F$'\t' -v OFS=$'\t' '$3!="." && $3!="" {gsub(/\..*$/,"",$1); gsub(/\..*$/,"",$3); print "'$sp'|"$3"\t'$sp'|"$1}' >> fish17.protein_to_gene.map
    else
        cat $gtp | awk -F$'\t' -v OFS=$'\t' '{gsub(/\.[0-9]*$/,"",$1); gsub(/\.[0-9]*$/,"",$2); print "'$sp'|"$2"\t'$sp'|"$1}' >> fish17.rna_to_gene.map
        cat $gtp | awk -F$'\t' -v OFS=$'\t' '$3!="." && $3!="" {gsub(/\.[0-9]*$/,"",$1); gsub(/\.[0-9]*$/,"",$3); print "'$sp'|"$3"\t'$sp'|"$1}' >> fish17.protein_to_gene.map
    fi
done

makeblastdb -in fish17.pep.fa  -dbtype prot
makeblastdb -in fish17.cdna.fa -dbtype nucl

mkdir -p fish17.blastp/
cmd=fish17.blastp/r01.blastp
file=$cmd.sh
>$file
echo '#!/bin/sh
#SBATCH --mem=8g -c 32 --time=240:00:00
module load blast/2.6.0+
np=`expr 2 \* $SLURM_CPUS_PER_TASK`
if ! [ -f all.m6.gz ]; then blastp -query ../fish17.pep.fa -db ../fish17.pep.fa -evalue 0.00001 -num_threads $np -outfmt "6 std btop qlen slen score gaps positive" -max_target_seqs 1000 | gzip -c > all.m6.gz; fi' >> $file
# cmd=r01.blastp; sbatch -J $cmd -o $cmd.o%A $cmd.sh

mkdir -p fish17.blastn/
cmd=fish17.blastn/r01.blastn
file=$cmd.sh
>$file
echo '#!/bin/sh
#SBATCH --mem=8g -c 32 --time=240:00:00
module load blast/2.6.0+
np=`expr 2 \* $SLURM_CPUS_PER_TASK`
if ! [ -f all.m6.gz ]; then blastn -task blastn -query ../fish17.cdna.fa -db ../fish17.cdna.fa -evalue 0.00001 -num_threads $np -outfmt "6 std btop qlen slen score gaps positive" -max_target_seqs 1000 | gzip -c > all.m6.gz; fi' >> $file
# cmd=r01.blastn; sbatch -J $cmd -o $cmd.o%A $cmd.sh

mybin=~/my_program3/src
for dir in fish17.blastn fish17.blastp
do
    cmd=$dir/r02.split_by_species
    cmd_file=$cmd.sh
    echo '#!/bin/sh' >$cmd_file
    if ! [ -d $dir/pairs ]; then mkdir $dir/pairs; fi
echo 'zcat all.m6.gz | perl -e '"'"'while(<>) { @t=split "\t",$_,-1; 
    $q=$t[0]; $t=$t[1]; 
    if ($q=~m/^([^\|]+)\|/) { $qsp=$1; }
    if ($t=~m/^([^\|]+)\|/) { $tsp=$1; }
    my $sp1 = $qsp;
    my $sp2 = $tsp;
    if (!exists $fh{$sp1}{$sp2}) { open $fh{$sp1}{$sp2}, "| gzip -c > pairs/$sp1.$sp2.m6.gz"; }
    print {$fh{$sp1}{$sp2}} join("\t",@t);
}
foreach my $sp1 (keys(%fh)) {
    foreach my $sp2 (keys(%{$fh{$sp1}})) {
        close $fh{$sp1}{$sp2};
    }
}'"'" >> $cmd_file;
done
# cmd=r02.split_by_species; sbatch --mem=8g --time=24:00:00 -J $cmd -o $cmd.o%A $cmd.sh

mybin=~/my_program3/src
for dir in fish17.blastn fish17.blastp
do
    if ! [ -d $dir/pairs.f ]; then mkdir $dir/pairs.f; fi
    if ! [ -d $dir/pairs.f2 ]; then mkdir $dir/pairs.f2; fi
    cmd=$dir/r03.join_colinear_simple
    cmd_file=$cmd.sh
    >$cmd_file
    if [ "$dir"=='fish17.blastn' ]; then iden=70; else iden=30; fi
    IFS=$'\n'
    for line1 in `cat fish.txt`
    do
        sp1=`echo $line1 | cut -f 4`
    for line2 in `cat fish.txt`
    do
        sp2=`echo $line2 | cut -f 4`
        n1=`echo $line1 | cut -f 9`;
        n2=`echo $line2 | cut -f 9`;
        if [ $n1 -gt 1 ] && [ $n2 -gt 1 ]
        then
            if [ $n1 -gt $n2 ]
            then
                n1=`expr $n1 / $n2`; n2=1;
            else
                n2=`expr $n2 / $n1`; n1=1;
            fi
        fi
echo 'perl '$mybin'/ohnolog/czl_blast_pair_filter.pl -i pairs/'$sp1'.'$sp2'.m6.gz -op pairs.f/'$sp1'.'$sp2'.f --piden-min '$iden' --ppositive-min 50 --e-max 0.001 --cov-min-min 50 --cov-max-min 75 --bit-frac-min 0.5 --piden-frac-min 0.75 --piden-d 5 --bit-frac-min2 0.3 --cov-min-min2 30 --cov-max-min2 50 --piden-d2 10 --copy '$sp1','$n1':'$sp2','$n2 >> $cmd_file
    done
    done
done
# cmd=r03.join_colinear_simple; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m blast,bedtools -b 20 -g 4 --partition quick --time=00:10:00 -J $cmd --logdir=$cmd.log -f $cmd.sh

for dir in fish17.blastn fish17.blastp
do
    if ! [ -d $dir/pairs.gene.f ]; then mkdir $dir/pairs.gene.f; fi
#    if ! [ -d $dir/pairs.gene.f2 ]; then mkdir $dir/pairs.gene.f2; fi
    cmd=$dir/r04.map_to_gene
    cmd_file=$cmd.sh
    >$cmd_file
    for sp1 in `cat fish.txt | cut -f 4`
    do
    for sp2 in `cat fish.txt | cut -f 4`
    do
if [ "$dir" == "fish17.blastn" ]; then map=../fish17.rna_to_gene.map; else map=../fish17.protein_to_gene.map; fi;
echo $mybin'/utility/czl_blast_map_id.pl -i pairs.f/'$sp1'.'$sp2'.f.join.m6.gz -o pairs.gene.f/'$sp1'.'$sp2'.f.join.m6.gz -m1 '$map' -m2 '$map' --remove-version' >> $cmd_file
#echo $mybin'/utility/czl_blast_map_id.pl -i pairs.f2/'$sp1'.'$sp2'.f2.join.m6.gz -o pairs.gene.f2/'$sp1'.'$sp2'.f2.join.m6.gz -m1 '$map' -m2 '$map' --remove-version' >> $cmd_file
    done
    done
done
# cmd=r04.map_to_gene; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m blast,bedtools -b 20 -g 4 --partition quick --time=00:10:00 -J $cmd --logdir=$cmd.log -f $cmd.sh

cat ../../carAur01/carAur01.gene.masked.gids  | awk '{print "GF|"$0}' > carAur01.gene.masked.sp_gids
cat ../carAur01/detect_het/run4/carAur01.noM.gene.masked.gids | awk '{print "GF|"$0}' > carAur01.gene.masked4.sp_gids
# remove masked gene
# {{{
for dir in fish17.blastn fish17.blastp
do
    if ! [ -d $dir/pairs.gene.f3 ]; then mkdir $dir/pairs.gene.f3; fi
    cmd=$dir/r05.remove_goldfish_het
    cmd_file=$cmd.sh
    >$cmd_file
    for sp1 in `cat fish.txt | cut -f 4`
    do
    for sp2 in `cat fish.txt | cut -f 4`
    do
if [ "$sp1" == "GF" ]
then
    if [ "$sp2" == "GF" ]
    then
echo 'sp1='$sp1'; sp2='$sp2'; zcat pairs.gene.f/$sp1.$sp2.f.join.m6.gz | sort -k1,1 -s > pairs.gene.f/$sp1.$sp2.f.join.m6.tmp;    join -t$'"'"'\t'"'"' -v 1 pairs.gene.f/$sp1.$sp2.f.join.m6.tmp ../carAur01.gene.masked.sp_gids | sort -k2,2 -s > pairs.gene.f/$sp1.$sp2.f.join.m6.tmp1;    join -t$'"'"'\t'"'"' -v 1 -1 2 -2 1 pairs.gene.f/$sp1.$sp2.f.join.m6.tmp1 ../carAur01.gene.masked.sp_gids | sort -k1,1 -s | gzip -c > pairs.gene.f3/$sp1.$sp2.f3.join.m6.gz;    rm pairs.gene.f/$sp1.$sp2.f.join.m6.tmp*;' >> $cmd_file
    else
echo 'sp1='$sp1'; sp2='$sp2'; zcat pairs.gene.f/$sp1.$sp2.f.join.m6.gz | sort -k1,1 -s > pairs.gene.f/$sp1.$sp2.f.join.m6.tmp;    join -t$'"'"'\t'"'"' -v 1 pairs.gene.f/$sp1.$sp2.f.join.m6.tmp ../carAur01.gene.masked.sp_gids | sort -k1,1 -s | gzip -c > pairs.gene.f3/$sp1.$sp2.f3.join.m6.gz;    rm pairs.gene.f/$sp1.$sp2.f.join.m6.tmp*;' >> $cmd_file
    fi
else
    if [ "$sp2" == "GF" ]
    then
echo 'sp1='$sp1'; sp2='$sp2'; zcat pairs.gene.f/$sp1.$sp2.f.join.m6.gz | sort -k2,2 -s | awk '"'"'{print $2"\t"$0}'"'"' > pairs.gene.f/$sp1.$sp2.f.join.m6.tmp;    join -t$'"'"'\t'"'"' -1 1 -2 1 -v 1 pairs.gene.f/$sp1.$sp2.f.join.m6.tmp ../carAur01.gene.masked.sp_gids | cut -f 2- | sort -k1,1 -s | gzip -c > pairs.gene.f3/$sp1.$sp2.f3.join.m6.gz;    rm pairs.gene.f/$sp1.$sp2.f.join.m6.tmp*;' >> $cmd_file
    else
echo 'sp1='$sp1'; sp2='$sp2'; ln -s '`cd $dir; pwd`'/pairs.gene.f/$sp1.$sp2.f.join.m6.gz pairs.gene.f3/$sp1.$sp2.f3.join.m6.gz;' >> $cmd_file
    fi
fi
    done
    done
done
# }}}

dir=fish17.blastn; cd $dir; zcat `ls pairs.gene.f3/??.??.f3.join.m6.gz | grep '\(ZF\|CC\|GF\|GC\|CF\).\(ZF\|CC\|GF\|GC\|CF\)'` | awk -F$'\t' -v OFS=$'\t' '{print $1,$2,"+",$12,$3}' | sort -k1,1 -k4,4nr > in5.edge; cd ../
