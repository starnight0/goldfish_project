module load R mafft

>fish4.gtp
>fish4.pep.fa
>fish4.cdna.fa
IFS=$'\n'
for line in `cat fish.txt | grep '\(ZF\|GC\|CC\|GF\)'`
do
    sp=`echo $line | cut -f 4`
    pep=`echo $line |  cut -f 5`
    cdna=`echo $line | cut -f 6`
    gtp=`echo $line | cut -f 8`
    zcat $pep  | sed -e '/^>/ s/\s.*$//' -e '/^>/ s/\..*$//' -e '/^>/ s/^>/>'$sp'|/' >> fish4.pep.fa
    zcat $cdna | sed -e '/^>/ s/\s.*$//' -e '/^>/ s/\..*$//' -e '/^>/ s/^>/>'$sp'|/' >> fish4.cdna.fa
    if [ "$sp" == "GC" ]
    then
        cat $gtp | awk -F$'\t' -v OFS=$'\t' '{gsub(/\..*$/,"",$1); $1="'$sp'|"$1; gsub(/\..*$/,"",$2); $2="'$sp'|"$2; gsub(/\..*$/,"",$3); $3="'$sp'|"$3; print $1,$2,$3}' >> fish4.gtp
    else
        cat $gtp | awk -F$'\t' -v OFS=$'\t' '{gsub(/\.[0-9]*$/,"",$1); $1="'$sp'|"$1; gsub(/\.[0-9]*$/,"",$2); $2="'$sp'|"$2; gsub(/\.[0-9]*$/,"",$3); $3="'$sp'|"$3; print $1,$2,$3}' >> fish4.gtp
    fi
done


mkdir -p fish4.cluster/run2
~/my_program3/src/ohnolog/phylogene.pl -a fish4.cluster/czl_ohno_syn.out3.4/anchor -gtp fish4.gtp -rna-fasta fish4.cdna.fa -pep-fasta fish4.pep.fa -chr-group fish4.cluster/chr_config2 -chr-group2 fish4.cluster/chr_config3 -c fish4.cluster/czl_ohno_syn.out3.4/cluster.txt -o fish4.cluster/run2/

#
mkdir gene_lost2
cat czl_ohno_syn.out3.4/rescue_m1.6.cluster.txt | perl -ne '
chomp; my @t=split /\t/; my @t1=split /,/,$t[2]; if (@t1<=1) {next;} my @c=(0,0,0,0); 
for my $a (@t1) {
    if ($a=~m/^ZF:([0-9]+)/) {$c[0]=$1;}
    elsif ($a=~m/^GC:([0-9]+)/) {$c[1]=$1;}
    elsif ($a=~m/^CC:([0-9]+)/) {$c[2]=$1;}
    elsif ($a=~m/^GF:([0-9]+)/) {$c[3]=$1;}
}
my $c = join "",@c;
if (!exists $fh{$c}) {
    open $fh{$c}, ">gene_lost2/rescue_m1.6.$c.id_to_clust.txt";
}
foreach my $aa (split /[,;]/,$t[3]) { if ($aa=~m/^\((.*)\)$/) { $aa=$1; $aa=~s/\|.*$//; } my @out=split /:/,$aa,4; my $out=join("\t",@out); $out=~s/:/\t/g; print {$fh{$c}} "$out\t$t[0]\t$t[1]\n";}
END { foreach my $c (keys(%fh)) { close $fh{$c}; } }
'
# ZF 
cat gene_lost2/rescue_m1.6.1???.id_to_clust.txt > gene_lost2/rescue_m1.6.1NNN.id_to_clust.txt
cat gene_lost2/rescue_m1.6.11??.id_to_clust.txt > gene_lost2/rescue_m1.6.11NN.id_to_clust.txt
cat gene_lost2/rescue_m1.6.10??.id_to_clust.txt > gene_lost2/rescue_m1.6.10NN.id_to_clust.txt
sp=ZF;
for cf in 1NNN 10NN 11NN
do
    f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt;
    if ! [ -d "gene_lost2/$cf.$sp" ]; then mkdir gene_lost2/$cf.$sp; fi
    cat $f1 | awk '$1~/'$sp'/' | cut -f 3 | sort -k1,1 > a1
    join -t$'\t' -j1 a1 /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/$sp/gene.ips.tsv > gene_lost2/$cf.$sp/gene.ips.tsv
    perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.$sp/gene.ips.tsv -o gene_lost2/$cf.$sp/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Interpro.txt
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Pfam.txt
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Gene3D.txt
done

# GC
cat gene_lost2/rescue_m1.6.?1??.id_to_clust.txt > gene_lost2/rescue_m1.6.N1NN.id_to_clust.txt
cat gene_lost2/rescue_m1.6.01??.id_to_clust.txt > gene_lost2/rescue_m1.6.01NN.id_to_clust.txt
sp=GC;
for cf in N1NN 01NN 11NN
do
f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt;
if ! [ -d "gene_lost2/$cf.$sp" ]; then mkdir gene_lost2/$cf.$sp; fi
cat $f1 | awk '$1~/'$sp'/' | cut -f 3 | sort -k1,1 > a1
join -t$'\t' -j1 a1 /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/$sp/gene.ips.tsv > gene_lost2/$cf.$sp/gene.ips.tsv
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.$sp/gene.ips.tsv -o gene_lost2/$cf.$sp/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Interpro.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Pfam.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Gene3D.txt
done

# compare ZF GC
mkdir -p gene_lost2/compare1/ZF.GC
mkdir -p gene_lost2/compare1/ZF2.ZF1
mkdir -p gene_lost2/compare1/GC2.GC1
for gs in Reactome KEGG Interpro Pfam GO2_CC.lv2 GO2_BP.lv2 GO2_MF.lv2 GO2_CC.lv3 GO2_BP.lv3 GO2_MF.lv3 GO2_CC.lv4 GO2_BP.lv4 GO2_MF.lv4 GO2_CC.lv5 GO2_BP.lv5 GO2_MF.lv5 GO2_CC.lv6 GO2_BP.lv6 GO2_MF.lv6
do
    echo $gs
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/1NNN.ZF/id_to_$gs.txt gene_lost2/N1NN.GC/id_to_$gs.txt gene_lost2/compare1/ZF.GC/$gs.enrich.txt
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/11NN.ZF/id_to_$gs.txt gene_lost2/10NN.ZF/id_to_$gs.txt gene_lost2/compare1/ZF2.ZF1/$gs.enrich.txt
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/11NN.GC/id_to_$gs.txt gene_lost2/01NN.GC/id_to_$gs.txt gene_lost2/compare1/GC2.GC1/$gs.enrich.txt
done
for gs in GO2_CC GO2_BP GO2_MF
do
>gene_lost2/compare1/GC2.GC1/$gs.enrich.txt
>gene_lost2/compare1/ZF2.ZF1/$gs.enrich.txt
for i in `seq 2 6`
do
cat gene_lost2/compare1/GC2.GC1/$gs.lv$i.enrich.txt | awk '{print $0"\t"'$i'}' >> gene_lost2/compare1/GC2.GC1/$gs.enrich.txt
cat gene_lost2/compare1/ZF2.ZF1/$gs.lv$i.enrich.txt | awk '{print $0"\t"'$i'}' >> gene_lost2/compare1/ZF2.ZF1/$gs.enrich.txt
done
done

# compare GF2 and GF12
cat gene_lost2/rescue_m1.6.[12][12]?2.id_to_clust.txt > gene_lost2/rescue_m1.6.BBN2.id_to_clust.txt
cat gene_lost2/rescue_m1.6.0[12]?2.id_to_clust.txt >> gene_lost2/rescue_m1.6.BBN2.id_to_clust.txt
cat gene_lost2/rescue_m1.6.[12]0?2.id_to_clust.txt >> gene_lost2/rescue_m1.6.BBN2.id_to_clust.txt
sp=GF;
cf=BBN2;
f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt;
mkdir gene_lost2/$cf.$sp
cat $f1 | awk '$1~/'$sp'/' | awk -F$'\t' 'n!=$5 {print; n=$5;}'| cut -f 3 | sort -k1,1 > a1
join -t$'\t' -j1 a1 /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/$sp/gene.ips.tsv > gene_lost2/$cf.$sp/gene.ips.tsv
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.$sp/gene.ips.tsv -o gene_lost2/$cf.$sp/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Interpro.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Pfam.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Gene3D.txt

cat gene_lost2/rescue_m1.6.[12][12]?[01].id_to_clust.txt > gene_lost2/rescue_m1.6.BBNA.id_to_clust.txt
cat gene_lost2/rescue_m1.6.0[12]?[01].id_to_clust.txt >> gene_lost2/rescue_m1.6.BBNA.id_to_clust.txt
cat gene_lost2/rescue_m1.6.[12]0?[01].id_to_clust.txt >> gene_lost2/rescue_m1.6.BBNA.id_to_clust.txt
sp=GF;
cf=BBNA;
f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt;
mkdir gene_lost2/$cf.$sp
cat $f1 | awk '$1~/'$sp'/' | cut -f 3 | sort -k1,1 > a1
join -t$'\t' -j1 a1 /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/$sp/gene.ips.tsv > gene_lost2/$cf.$sp/gene.ips.tsv
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.$sp/gene.ips.tsv -o gene_lost2/$cf.$sp/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Interpro.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Pfam.txt

cat gene_lost2/rescue_m1.6.[12][12]??.id_to_clust.txt > gene_lost2/rescue_m1.6.BBNN.id_to_clust.txt
cat gene_lost2/rescue_m1.6.0[12]??.id_to_clust.txt >> gene_lost2/rescue_m1.6.BBNN.id_to_clust.txt
cat gene_lost2/rescue_m1.6.[12]0??.id_to_clust.txt >> gene_lost2/rescue_m1.6.BBNN.id_to_clust.txt
sp=GF;
cf=BBNN;
f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt;
mkdir gene_lost2/$cf.$sp
cat $f1 | awk '$1~/'$sp'/' | awk -F$'\t' 'n!=$5 {print; n=$5;}'| cut -f 3 | sort -k1,1 > a1
join -t$'\t' -j1 a1 /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/$sp/gene.ips.tsv > gene_lost2/$cf.$sp/gene.ips.tsv
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.$sp/gene.ips.tsv -o gene_lost2/$cf.$sp/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Interpro.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Pfam.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Gene3D.txt

cat gene_lost2/rescue_m1.6.[12][12]?1.id_to_clust.txt > gene_lost2/rescue_m1.6.BBN1.id_to_clust.txt
cat gene_lost2/rescue_m1.6.0[12]?1.id_to_clust.txt >> gene_lost2/rescue_m1.6.BBN1.id_to_clust.txt
cat gene_lost2/rescue_m1.6.[12]0?1.id_to_clust.txt >> gene_lost2/rescue_m1.6.BBN1.id_to_clust.txt
sp=GF;
cf=BBN1;
f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt;
mkdir gene_lost2/$cf.$sp
cat $f1 | awk '$1~/'$sp'/' | cut -f 3 | sort -k1,1 > a1
join -t$'\t' -j1 a1 /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/$sp/gene.ips.tsv > gene_lost2/$cf.$sp/gene.ips.tsv
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.$sp/gene.ips.tsv -o gene_lost2/$cf.$sp/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Interpro.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Pfam.txt
cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Gene3D.txt
# compare GF2 and GF12
mkdir -p gene_lost2/compare1/GF.GF.1
mkdir -p gene_lost2/compare1/GF.GF.2
for gs in Reactome KEGG Interpro Pfam GO2_CC.lv2 GO2_BP.lv2 GO2_MF.lv2 GO2_CC.lv3 GO2_BP.lv3 GO2_MF.lv3 GO2_CC.lv4 GO2_BP.lv4 GO2_MF.lv4 GO2_CC.lv5 GO2_BP.lv5 GO2_MF.lv5 GO2_CC.lv6 GO2_BP.lv6 GO2_MF.lv6
do
    echo $gs
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/BBN2.GF/id_to_$gs.txt gene_lost2/BBNN.GF/id_to_$gs.txt gene_lost2/compare1/GF.GF.1/$gs.enrich.txt
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/BBN2.GF/id_to_$gs.txt gene_lost2/BBN1.GF/id_to_$gs.txt gene_lost2/compare1/GF.GF.2/$gs.enrich.txt
done

# GF lost, use ZF as ref
cat gene_lost2/rescue_m1.6.1??0.id_to_clust.txt > gene_lost2/rescue_m1.6.1NN0.id_to_clust.txt
cat gene_lost2/rescue_m1.6.1??1.id_to_clust.txt > gene_lost2/rescue_m1.6.1NN1.id_to_clust.txt
cat gene_lost2/rescue_m1.6.1??2.id_to_clust.txt > gene_lost2/rescue_m1.6.1NN2.id_to_clust.txt
sp=ZF;
for cf in 1NN0 1NN1 1NN2
do
    f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt;
    mkdir gene_lost2/$cf.$sp
    cat $f1 | awk '$1~/'$sp'/' | cut -f 3 | sort -k1,1 > a1
    join -t$'\t' -j1 a1 /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/$sp/gene.ips.tsv > gene_lost2/$cf.$sp/gene.ips.tsv
    perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.$sp/gene.ips.tsv -o gene_lost2/$cf.$sp/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Interpro.txt
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Pfam.txt
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Gene3D.txt
done

mkdir gene_lost2/compare1/GF.GF.3
mkdir gene_lost2/compare1/GF.GF.4
for gs in Reactome KEGG Interpro Pfam GO2_CC.lv2 GO2_BP.lv2 GO2_MF.lv2 GO2_CC.lv3 GO2_BP.lv3 GO2_MF.lv3 GO2_CC.lv4 GO2_BP.lv4 GO2_MF.lv4 GO2_CC.lv5 GO2_BP.lv5 GO2_MF.lv5 GO2_CC.lv6 GO2_BP.lv6 GO2_MF.lv6
do
    echo $gs
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/1NN2.ZF/id_to_$gs.txt gene_lost2/1NN1.ZF/id_to_$gs.txt gene_lost2/compare1/GF.GF.3/$gs.enrich.txt
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/1NN2.ZF/id_to_$gs.txt gene_lost2/1NN0.ZF/id_to_$gs.txt gene_lost2/compare1/GF.GF.4/$gs.enrich.txt
done
# GF use GC as reference
cat gene_lost2/rescue_m1.6.?1?0.id_to_clust.txt > gene_lost2/rescue_m1.6.N1N0.id_to_clust.txt
cat gene_lost2/rescue_m1.6.?1?1.id_to_clust.txt > gene_lost2/rescue_m1.6.N1N1.id_to_clust.txt
cat gene_lost2/rescue_m1.6.?1?2.id_to_clust.txt > gene_lost2/rescue_m1.6.N1N2.id_to_clust.txt
sp=GC;
for cf in N1N0 N1N1 N1N2
do
    f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt;
    mkdir gene_lost2/$cf.$sp
    cat $f1 | awk '$1~/'$sp'/' | cut -f 3 | sort -k1,1 > a1
    join -t$'\t' -j1 a1 /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/big/ips/$sp/gene.ips.tsv > gene_lost2/$cf.$sp/gene.ips.tsv
    perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.$sp/gene.ips.tsv -o gene_lost2/$cf.$sp/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Interpro.txt
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Pfam.txt
    cat gene_lost2/$cf.$sp/gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.$sp/id_to_Gene3D.txt
done

mkdir gene_lost2/compare1/GC.GF2.GF1
mkdir gene_lost2/compare1/GC.GF2.GF0
for gs in Reactome KEGG Interpro Pfam GO2_CC.lv2 GO2_BP.lv2 GO2_MF.lv2 GO2_CC.lv3 GO2_BP.lv3 GO2_MF.lv3 GO2_CC.lv4 GO2_BP.lv4 GO2_MF.lv4 GO2_CC.lv5 GO2_BP.lv5 GO2_MF.lv5 GO2_CC.lv6 GO2_BP.lv6 GO2_MF.lv6
do
    echo $gs
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/N1N2.GC/id_to_$gs.txt gene_lost2/N1N1.GC/id_to_$gs.txt gene_lost2/compare1/GC.GF2.GF1/$gs.enrich.txt
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/N1N2.GC/id_to_$gs.txt gene_lost2/N1N0.GC/id_to_$gs.txt gene_lost2/compare1/GC.GF2.GF0/$gs.enrich.txt
done

for gs in GO2_CC GO2_BP GO2_MF
do
>gene_lost2/compare1/GC.GF2.GF1/$gs.enrich.txt
>gene_lost2/compare1/GC.GF2.GF0/$gs.enrich.txt
for i in `seq 2 6`
do
cat gene_lost2/compare1/GC.GF2.GF1/$gs.lv$i.enrich.txt | awk '{print $0"\t"'$i'}' >> gene_lost2/compare1/GC.GF2.GF1/$gs.enrich.txt
cat gene_lost2/compare1/GC.GF2.GF0/$gs.lv$i.enrich.txt | awk '{print $0"\t"'$i'}' >> gene_lost2/compare1/GC.GF2.GF0/$gs.enrich.txt
done
done

# merge all levels of GO
for dir in GC.GF2.GF0 GC.GF2.GF1 GF.GF.3 GF.GF.4
do
    for gs in GO2_CC GO2_BP GO2_MF
    do
        dir1=gene_lost2/compare1/$dir;
        cat $dir1/$gs.lv2.enrich.txt | head -n 1 | awk '{print $0"\tGO_level"}' > $dir1/$gs.enrich.txt
        for i in `seq 2 6`
        do
            cat $dir1/$gs.lv$i.enrich.txt | tail -n +2 | awk '{print $0"\t"'$i'}' >> $dir1/$gs.enrich.txt
        done
    done
    for gs in GO2_CC GO2_BP GO2_MF KEGG Reactome Interpro Pfam
    do
        dir1=gene_lost2/compare1/$dir;
        cat $dir1/$gs.enrich.txt | awk -F$'\t' 'NR==1 {print} NR>1 { n=0; for (i=0;i<1;i++) {j1=8+i*10; j2=9+i*10; if ($j1<0.01 && $j2<0.1) {n++;} } if (n>0) {print}}' > $dir1/$gs.enrich.f.txt
    done
done


#

for cf in 1112 1121 '11--'
do
    f1=gene_lost2/rescue_m1.6.$cf.id_to_clust.txt
    cat $f1 | awk '$1~/ZF/' | cut -f 3 | sort -k1,1 > a1
    join -t$'\t' -j1 a1 ~/data/ensembl85/ips/Danio_rerio.GRCz10.pep.all.ips.parse/gene.ips.tsv > gene_lost2/$cf.ZF.gene.ips.tsv
    if ! [ -d gene_lost2/$cf.ZF.gene.ips.parse ]; then mkdir gene_lost2/$cf.ZF.gene.ips.parse; fi
    perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene_lost2/$cf.ZF.gene.ips.tsv -o gene_lost2/$cf.ZF.gene.ips.parse/ -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
    cat gene_lost2/$cf.ZF.gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.ZF.gene.ips.parse/id_to_Interpro.txt
    cat gene_lost2/$cf.ZF.gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.ZF.gene.ips.parse/id_to_Pfam.txt
    cat gene_lost2/$cf.ZF.gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > gene_lost2/$cf.ZF.gene.ips.parse/id_to_Gene3D.txt
done
# compare
mkdir gene_lost2/compare;
for cf in 1112 1121
do
for d in GO2_MF GO2_CC GO2_BP
do
echo $cf $d
for i in 1 2 3 4
do
#    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/$cf.ZF.gene.ips.parse/id_to_$d.lv$i.txt ~/data/ensembl85/ips/Danio_rerio.GRCz10.pep.all.ips.parse/gene.ips.id_to_$d.lv$i.txt gene_lost2/compare/$cf.$d.lv$i.enrich.txt
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/$cf.ZF.gene.ips.parse/id_to_$d.lv$i.txt gene_lost2/11--.ZF.gene.ips.parse/id_to_$d.lv$i.txt gene_lost2/compare/$cf.$d.lv$i.enrich.txt
done
done
done

for cf in 1112 1121
do
for d in Reactome KEGG Interpro Pfam
do
echo $cf $d
#    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/$cf.ZF.gene.ips.parse/id_to_$d.txt ~/data/ensembl85/ips/Danio_rerio.GRCz10.pep.all.ips.parse/gene.ips.id_to_$d.txt gene_lost2/compare/$cf.$d.enrich.txt
    Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript gene_lost2/$cf.ZF.gene.ips.parse/id_to_$d.txt gene_lost2/11--.ZF.gene.ips.parse/id_to_$d.txt gene_lost2/compare/$cf.$d.enrich.txt
done
done
