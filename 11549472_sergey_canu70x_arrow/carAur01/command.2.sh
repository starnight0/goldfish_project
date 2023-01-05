# when move files to datashare, DO:
# setfacl --modify user:webcpu:r-- FILE
###########################################

exit 0;


kent=/data/genome/jksrc_v352/kent/src/hg/lib/
asm1_dir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01
asm1=carAur01

mkdir big/est2genome
mkdir bbi/est2genome


##### copy readToTig #########
liftOver ../../asm.contigs.readToTig.bed ../asm_to_arrow.liftOver.chain big/carAur01.contigs.readToTig.bed  big/carAur01.contigs.readToTig.unmap 
cat big/carAur01.contigs.readToTig.bed | grep -v tig00037772 > big/carAur01.contigs.readToTig.bed
bedToBigBed carAur01.contigs.readToTig.bed carAur01.chromSizes bbi/carAur01.contigs.readToTig.bb
bedToBigBed big/carAur01.noM.contigs.readToTig.bed carAur01.noM.sm.fa.fai bbi/carAur01.noM.contigs.readToTig.bb
bedtools genomecov -i carAur01.contigs.readToTig.bed -g carAur01.chromSizes -bg | gzip -c > carAur01.contigs.readToTig.bg.gz
bedGraphToBigWig carAur01.contigs.readToTig.bg carAur01.chromSizes carAur01.contigs.readToTig.bw
cat big/carAur01.contigs.readToTig.bg | grep -v tig00037772 > big/carAur01.noM.contigs.readToTig.bg 
bedGraphToBigWig big/carAur01.noM.contigs.readToTig.bg carAur01.noM.sm.fa.fai bbi/carAur01.noM.contigs.readToTig.bw
# trim end 500bp
cat big/carAur01.noM.contigs.readToTig.bed | awk -F$'\t' -v OFS=$'\t' '$3-$2-1000>0{$2+=500; $3-=500; print $0}' > big/carAur01.noM.contigs.readToTig.trim_end_500bp.bed 
bedtools genomecov -i big/carAur01.noM.contigs.readToTig.trim_end_500bp.bed -g carAur01.withM.sm.fa.fai -bg > big/carAur01.noM.contigs.readToTig.trim_end_500bp.bg
bedGraphToBigWig big/carAur01.noM.contigs.readToTig.trim_end_500bp.bg carAur01.noM.sm.fa.fai bbi/carAur01.noM.contigs.readToTig.trim_end_500bp.bw
# trim end 1000bp
cat big/carAur01.noM.contigs.readToTig.bed | awk -F$'\t' -v OFS=$'\t' '$3-$2>=3000 {$2+=1000; $3-=1000; print $0}' > big/carAur01.noM.contigs.readToTig.trim_end_1000bp.bed 
bedtools genomecov -i big/carAur01.noM.contigs.readToTig.trim_end_1000bp.bed -g carAur01.withM.sm.fa.fai -bg > big/carAur01.noM.contigs.readToTig.trim_end_1000bp.bg
bedGraphToBigWig big/carAur01.noM.contigs.readToTig.trim_end_1000bp.bg carAur01.noM.sm.fa.fai bbi/carAur01.noM.contigs.readToTig.trim_end_1000bp.bw

##### short-read map #######
cp ../HC2NCBCXX_bwa/bedtools.genomecov.Q0.bw bbi/HC2NCBCXX_bwa.bedtools.genomecov.Q0.bw 
cp ../HC2NCBCXX_bwa/bedtools.genomecov.Q3.bw bbi/HC2NCBCXX_bwa.bedtools.genomecov.Q3.bw 
cp ../HC2NCBCXX_bwa/PP.bedtools.genomecov.Q3.bw bbi/HC2NCBCXX_bwa.PP.bedtools.genomecov.Q3.bw 

##### copy tigInfo #########
cat ../../asm.contigs.layout.tigInfo | awk -F$'\t' -v OFS=$'\t' '$6=="contig" {$1=sprintf("tig%08i_arrow", $1); print}' > carAur01.ctgInfo.tmp
head -n 1 ../../asm.contigs.layout.tigInfo  >  carAur01.ctgInfo.header.txt
join -t$'\t' -j 1 carAur01.chromSizes carAur01.ctgInfo.tmp | cut -f 1,2,4- > carAur01.ctgInfo
rm  carAur01.ctgInfo.tmp
#dp=`cat carAur01.ctgInfo  | awk -v n=0 -v l=0 '$5>=20 && $5<=100 {l+=$2; n+=$2*$5} END{print n/l}'`
dp=45
cat carAur01.ctgInfo  | awk -F$'\t' -v dp=$dp '{print $0"\t"$5/dp}' > carAur01.ctgInfo2
cat carAur01.ctgInfo2 | awk -F$'\t' '$10>=0.6 && $10<1.6 {print $1}' > carAur01.norm_coverage.ctg_id
cat carAur01.ctgInfo2 | awk -F$'\t' -v OFS=$'\t' -v n0=0 -v n1=0 -v n2=0 '{if ($10<0.6) {n0++; m0+=$2;} else if ($10<1.6) {n1++;m1+=$2;} else {n2++;m2+=$2;}} END{print "Coverage","Half","One","Multiple"; print "N_CTG",n0,n1,n2; print "total_bp",m0,m1,m2}' > carAur01.ctgInfo2.stat

echo >>  carAur01.ctgInfo2.stat
cat carAur01.ctgInfo2 | awk -F$'\t' -v OFS=$'\t' -v n0=0 -v n1=0 '{if ($7=="no") {n0++; m0+=$2;} else {n1++;m1+=$2;} } END{print "Type","Norm","Repeat"; print "N_CTG",n0,n1; print "total_bp",m0,m1}' >> carAur01.ctgInfo2.stat

echo >>  carAur01.ctgInfo2.stat
echo "NL50 of Half-coverage CTG" >> carAur01.ctgInfo2.stat
cat carAur01.ctgInfo2  | awk -F$'\t' '$10<0.6 {print $1"\t"$2}' | ~/my_program3/src/assembler/NL50.pl - >>  carAur01.ctgInfo2.stat

echo >>  carAur01.ctgInfo2.stat
echo "NL50 of Normal-coverage CTG" >> carAur01.ctgInfo2.stat
cat carAur01.ctgInfo2  | awk -F$'\t' '$10>=0.6 && $10<1.6 {print $1"\t"$2}' | ~/my_program3/src/assembler/NL50.pl - >>  carAur01.ctgInfo2.stat
############################

##### combine contig stats ##########
~/my_program3/src/utility/czl_tab_join.pl -i1 contig.stat -i2 carAur01.not_contained.fasta.fai -1 2 -2 1 -o - | awk '$10!="."' > contig.not_contained.stat
sort -k2,2 contig.stat > contig.sorted.stat
join -t$'\t' -1 2 -2 1 contig.sorted.stat ../asm.contigs.layout.ctgInfo2 | cut -f 1,3- > contig.stat2
join -t$'\t' -1 1 -2 1 contig.contained.stat ../asm.contigs.layout.ctgInfo2 | cut -f 1,3- > contig.contained.stat2
##################

#######################################################
#  gene track
#######################################################
# {{{
#********************************************
# add 2018-03-16
#********************************************
cp ../maker4b/carAur01.add_more_name_ips.gff ./big/carAur01.gene.gff
gff3ToGenePred -defaultCdsStatusToUnknown -attrsOut=big/carAur01.gene.gff.attr big/carAur01.gene.gff big/carAur01.gene.gp
#
cat big/carAur01.gene.gff | perl -ne '
if (m/^#/ || m/^\s*$/) { next; } chomp; my @t=split "\t",$_,-1;
my %info; foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa,2; $info{$u}=$v;}
if ($t[2]=~m/gene/) {$gname=$info{Name};}
elsif ($t[2]=~m/(RNA|transcript)/) { print "$info{ID}\t$info{Name}\t$gname\n"; } ' > big/carAur01.gene.gff.names
#
cat big/carAur01.gene.gff | perl -ne '
if (m/^#/ || m/^\s*$/) { next; } chomp; my @t=split "\t",$_,-1;
my %info; foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa,2; $info{$u}=$v;}
if ($t[2]=~m/(RNA|transcript)/) {
    $out[0]=$info{ID};
    if ($info{cds_stat}=~m/5/) { $out[1]="cmpl"; } else { $out[1]="incmpl"; }
    if ($info{cds_stat}=~m/3/) { $out[2]="cmpl"; } else { $out[2]="incmpl"; }
    print join("\t",@out[0..2]), "\n";
} ' | sort -k1,1 > big/carAur01.gene.gff.cds_stat
#cat big/carAur01.gene.gp | cut -f 1,15 | sort -k1,1 > a2
#join -t$'\t' -j1 a1 a2 > big/carAur01.gene.gff.knowCds

genePredToBigGenePred -geneNames=big/carAur01.gene.gff.names big/carAur01.gene.gp stdout |  awk -F$'\t' '{print $4"\t"$0}' | sort -k1,1 > a
join -t$'\t' -1 1 -2 1 a big/carAur01.gene.gff.cds_stat | cut -f 2- | awk -F$'\t' -v OFS=$'\t' '{g=$13;$13=$18;$18=g;$14=$21;$15=$22;NF=20; print $0}' | sort -k1,1 -k2,2n > big/carAur01.gene.bgp
join -t$'\t' carAur01.norm_coverage.ctg_id big/carAur01.noM.gene.bgp > big/carAur01.noM.gene.norm_coverage.bgp
cat big/carAur01.noM.gene.norm_coverage.bgp | cut -f 18 | sort | uniq > big/carAur01.noM.gene.norm_coverage.bgp.gid
cat big/carAur01.gene.bgp | grep -v tig00037772 > big/carAur01.noM.gene.bgp
bedToBigBed -extraIndex=name,name2,geneName,geneName2 -tab -as=$kent/bigGenePred.as -type=bed12+8 big/carAur01.noM.gene.bgp carAur01.withM.sm.fa.fai bbi/carAur01.noM.gene.bb

cat big/carAur01.gene.gff | perl -ne '
if (m/^#/ || m/^\s*$/) { next; } chomp; my @t=split "\t",$_,-1;
my %info; foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa,2; $info{$u}=$v;}
if ($t[2]=~m/(gene)/) { print join("\t", ($t[0], $t[3]-1, $t[4], $info{ID}, $t[5], $t[6]) ), "\n" }' | sort -k1,1 -k2,2n > big/carAur01.gene.bed

cat big/carAur01.gene.bgp | cut -f  1-6 > big/carAur01.transcript.bed

cat big/carAur01.gene.gff | perl -ne '
if (m/^#/ || m/^\s*$/) { next; } chomp; my @t=split "\t",$_,-1;
my %info; foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa,2; $info{$u}=$v;}
if ($t[2]=~m/exon/) { print join("\t", ($t[0], $t[3]-1, $t[4], $info{ID}, $t[5], $t[6]) ), "\n" }' | sort -k1,1 -k2,2n > big/carAur01.exon.bed

cat big/carAur01.gene.gff | perl -ne '
if (m/^#/ || m/^\s*$/) { next; } chomp; my @t=split "\t",$_,-1;
my %info; foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa,2; $info{$u}=$v;}
if ($t[2]=~m/CDS/i) { print join("\t", ($t[0], $t[3]-1, $t[4], $info{ID}, $t[5], $t[6]) ), "\n" }' | sort -k1,1 -k2,2n > big/carAur01.cds.bed

#************ Fetch exon sequence ************
gffread big/carAur01.noM.gene.gff -g carAur01.withM.sm.fa --force-exons -w  big/carAur01.noM.exon.fa -x big/carAur01.noM.cds.fa -y big/carAur01.noM.tr_cds.fa -W 

#*****************************************
# add 2018-03-20
#*****************************************
cat big/carAur01.gene.gff | perl -ne '
if (m/^#/ || m/^\s*$/) { next; } chomp; my @t=split "\t",$_,-1;
my %info; foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa,2; $info{$u}=$v;}
if ($t[2]=~m/gene/) {
    if (exists $tr{tid}) { print "$tr{chr}\t$tr{begin}\t$tr{end}\t$tr{gid}|$tr{tid}|$tr{pid}|$tr{gname}|$tr{exon_len}|$tr{cds_len}|$tr{tname}|$tr{gbiotype}|$tr{tbiotype}|$tr{exon_num}|$tr{cds_num}\t1000\t$tr{strand}\n";  undef %tr; }
    $tr{gid}=$info{ID}; $tr{gname}=$info{Name}; $first=1;}
elsif ($t[2]=~/RNA/) {
if (exists $tr{tid}) { print "$tr{chr}\t$tr{begin}\t$tr{end}\t$tr{gid}|$tr{tid}|$tr{pid}|$tr{gname}|$tr{exon_len}|$tr{cds_len}|$tr{tname}|$tr{gbiotype}|$tr{tbiotype}|$tr{exon_num}|$tr{cds_num}\t1000\t$tr{strand}\n";  }
if ($first) { delete $tr{gbiotype}; $first=0; }
$tr{chr} = $t[0];  $tr{begin}=$t[3]-1;  $tr{end} = $t[4];  $tr{strand} = $t[6];  
$tr{tid}=$info{ID};  $tr{tname}=$info{Name};
if ($t[2]=~/mRNA/) { $tr{pid}=$tr{tid}; } else { $tr{pid}="."; }
$tr{tbiotype} = $t[2];  if (!exists $tr{gbiotype}) { $tr{gbiotype} = $t[2]; }
$tr{exon_num}=$tr{cds_num}=0;   $tr{exon_len}=$tr{cds_len}=0;
} elsif ($t[2]=~/exon/) {
$tr{exon_num}++;  $tr{exon_len}+=$t[4]-$t[3]+1;
} elsif ($t[2]=~/cds/i) { 
$tr{cds_num}++;   $tr{cds_len}+=$t[4]-$t[3]+1;
}       
END {
    if (exists $tr{tid}) { print "$tr{chr}\t$tr{begin}\t$tr{end}\t$tr{gid}|$tr{tid}|$tr{pid}|$tr{gname}|$tr{exon_len}|$tr{cds_len}|$tr{tname}|$tr{gbiotype}|$tr{tbiotype}|$tr{exon_num}|$tr{cds_num}\t1000\t$tr{strand}\n"; }
}' > big/carAur01.gene.gff.bed
cat big/carAur01.gene.gff.bed | cut -f 4 | awk -F'|' -v OFS=$'\t' '{print $1,$2,$3,$4,$7,$8,$9,$10,$11}' > big/carAur01.gene.gff.gtpnnbbcc
#*****************************************

#---------------------------------------
# gene a .cds file for pslToBigPsl 
#---------------------------------------
cat big/carAur01.noM.gene.bgp | perl -ne '
if (m/^#/ || m/^\s*$/) { next; }  chomp;
my @t=split "\t",$_;  my $b; my $e;
$cds_start=$t[6];   $cds_end=$t[7];
if ($t[9]==1) {
    $b = $t[1]+$t[11];    $e = $t[1]+$t[11]+$t[10];
    if ($b<$cds_start) { $b = $cds_start; }
    if ($e>$cds_end) { $e = $cds_end; }
    print $t[0] , "\t" , $b+1 , ".." , $e , "\n";
} else {
    my @sz = split ",", $t[10];
    my @b = split ",", $t[11];
    my @aa = ();
    for (my $i=0; $i<=$#sz; $i++) {
        $b = $t[1]+$b[$i];    $e = $b+$sz[$i];
        if ($e <= $cds_start) {next;}
        if ($b >= $cds_end) {next;}
        if ($b<$cds_start) { $b = $cds_start; }
        if ($e>$cds_end) { $e = $cds_end; }
        push @aa, $b+1 . ".." . $e; 
    }
    if ($#aa==0) { print "$t[0]\t$aa[0]\n"; }
    else { print "$t[0]\tjoin(", join(",", @aa), ")", "\n"; }
}' > big/carAur01.noM.CDS.for_bigpsl
#---------------------------------------
#}}}
#######################################################


##########################
# get a list of het genes
##########################
# {{{
join -t$'\t' -j 1 contig.contained.names big/carAur01.noM.gene.bgp > big/carAur01.noM.gene.masked1.bgp
cat big/carAur01.noM.gene.masked1.bgp | cut -f 4 | sort | uniq > carAur01.gene.masked1.tids
cat big/carAur01.noM.gene.masked1.bgp | cut -f 18| sort | uniq > carAur01.gene.masked1.gids
zcat ../maker4a/t.blast/sp5.pairs.2/carAur.carAur.f.join.m6.gz | sed 's/carAur|//g' > ./carAur.carAur.f.join.no_sp.blastn.m6
#~/my_program3/src/utility/czl_tab_join.pl -i1 carAur.carAur.f.join.no_sp.blastn.m6 -i2 carAur01.gene.annot.contained.tids -1 1 -2 1 -o - | awk -F$'\t' '$1!="." && $NF=="."' | cut -f 1-20 | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 carAur01.gene.annot.contained.tids -1 2 -2 1 -o - | awk -F$'\t' '$1!="." && $NF=="."' | cut -f 1-20  > carAur.carAur.f.join.no_sp.blastn.no_contained_ctg.m6
#~/my_program3/src/annot_genome/czl_goldfish_het_gene.pl --bgp carAur01.gene.annot.bgp --ctg-info ../carAur03/carAur01.ctgInfo3 --self-m6 carAur.carAur.f.join.no_sp.blastn.no_contained_ctg.m6 -o carAur01.
#cat contig.contained.names carAur01.het_remove_contig.txt | sort | uniq > carAur01.to_remove_contig.txt
cut -f 1 detect_het/run2/het_remove_gene.txt | sort | uniq > carAur01.gene.masked2.gids
cat carAur01.gene.masked1.gids carAur01.gene.masked2.gids | sort -k1,1 | uniq > carAur01.gene.masked.gids
cut -f 1,2 big/carAur01.gene.gff.gtpnnbbcc | sort -k1,1 > a
join -t$'\t' -j 1 a carAur01.gene.masked2.gids | cut -f 2 | sort > carAur01.gene.masked2.tids
join -t$'\t' -j 1 a carAur01.gene.masked.gids | cut -f 2 | sort > carAur01.gene.masked.tids

cat big/carAur01.noM.gene.bgp | awk -F$'\t' '{print $18"\t"$0}'| sort -k1,1 > a 
join -t$'\t' -j1 carAur01.gene.unmasked.gids a | cut -f 2-21 | sort -k1,1 -k2,2n > big/carAur01.noM.gene.unmasked.bgp
join -t$'\t' -j1 carAur01.gene.masked.gids a | cut -f 2-21 | sort -k1,1 -k2,2n > big/carAur01.noM.gene.masked.bgp

cut -f 1 big/carAur01.gene.gff.gtpnnbbcc | sort | uniq > carAur01.gene.gids
cut -f 2 big/carAur01.gene.gff.gtpnnbbcc | sort | uniq > carAur01.gene.tids
cat carAur01.gene.gids carAur01.gene.masked.gids | sort | uniq -u > carAur01.gene.unmasked.gids 
cat carAur01.gene.tids carAur01.gene.masked.tids | sort | uniq -u > carAur01.gene.unmasked.tids 

# produce a unmasked gene gff
cat big/carAur01.gene.gff | perl -e '
my %gid;
open IN, "<carAur01.gene.masked.gids";
while (<IN>) {if (/^#/ || /^\s*$/) {next;} chomp; my @t=split /\t/,$_; $gid{$t[0]}++;}
close IN;
open OUT1, ">big/carAur01.gene.masked.gff";
open OUT2, ">big/carAur01.gene.unmasked.gff";
my $filt=0;
while(<>) {
    if (/^#/ || /^\s*$/) {next;}
    chomp; my @t=split /\t/,$_;
    if ($t[2] eq "gene") {
        my %info;
        foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa; $info{$u}=$v;}
        if (exists $gid{$info{ID}}) { $filt=1; print OUT1 $_, "\n"; } else { $filt=0; print OUT2 $_, "\n";}
    } elsif (!$filt) {
        print OUT2 $_, "\n";
	} else {
        print OUT1 $_, "\n";
    }
} close OUT1;  close OUT2; '
cat big/carAur01.gene.unmasked.gff | grep -v tig00037772 > big/carAur01.noM.gene.unmasked.gff
cat big/carAur01.gene.masked.gff | grep -v tig00037772 > big/carAur01.noM.gene.masked.gff

# get a list of masked contigs
cat detect_het/het_remove_contig.txt | cut -f 1 > a1
echo tig00037772_arrow >> a1
cat contig.contained.names a1 | sort -k1,1 | uniq > carAur01.masked_ctg_id
cut -f 1 carAur01.withM.sm.fa.fai > a1
cat a1 carAur01.masked_ctg_id | sort -k1,1 | uniq -u | grep -v tig00037772 > carAur01.unmasked_ctg_id

cut -f 1,2 carAur01.withM.sm.fa.fai | sort -k1,1 > a2
join -t$'\t' -j1 carAur01.unmasked_ctg_id a2 > carAur01.unmasked_ctg_id_size; 
cat carAur01.unmasked_ctg_id_size | awk '{print $1"\t1\t"$2}' > carAur01.unmasked_ctg_id.bed3; 

join -t$'\t' -j 1 carAur01.masked_ctg_id big/carAur01.gene.bgp | cut -f 4 | sort | uniq > carAur01.gene.masked3.tids
join -t$'\t' -j 1 carAur01.masked_ctg_id big/carAur01.gene.bgp | cut -f 18 | sort | uniq > carAur01.gene.masked3.gids
cat carAur01.gene.masked[123].tids | sort | uniq >  carAur01.gene.masked123.tids
cat carAur01.gene.masked[123].gids | sort | uniq >  carAur01.gene.masked123.gids
~/my_program3/src/utility/czl_fasta_fetch -i carAur01.withM.sm.fa -o carAur02.withM.sm. -m 1 --list carAur01.unmasked_ctg_id
samtools faidx carAur02.withM.sm.fa 
cat big/carAur01.noM.gene.bgp | awk -F$'\t' '{print $18"\t"$0}'| sort -k1,1 > a 
join -t$'\t' -j1 -v 1 a carAur01.gene.masked123.gids | cut -f 2-21 | sort -k1,1 -k2,2n > big/carAur02.noM.gene.bgp
cat big/carAur02.noM.gene.bgp | cut -f 18 | sort | uniq > carAur02.noM.gene.gids
cat big/carAur02.noM.gene.bgp | cut -f 4  | sort | uniq > carAur02.noM.gene.tids
# produce a unmasked gene gff for carAur02
cat big/carAur01.gene.gff | perl -e '
my %gid;
open IN, "<carAur01.gene.masked123.gids";
while (<IN>) {if (/^#/ || /^\s*$/) {next;} chomp; my @t=split /\t/,$_; $gid{$t[0]}++;}
close IN;
my $filt=0;
while(<>) {
    if (/^#/ || /^\s*$/) {next;}
    chomp; my @t=split /\t/,$_;
    if ($t[2] eq "gene") {
        my %info;
        foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa; $info{$u}=$v;}
        if (exists $gid{$info{ID}}) { $filt=1; } else { $filt=0; print $_, "\n";}
    } elsif (!$filt) {
        print $_, "\n";
    }
}
' > big/carAur02.noM.gene.gff


# ../maker4a/p.blast/sp5.pairs.2/carAur.carAur.f2.join.no_sp.m6.gz 
# OUTPUT:  carAur01.het_remove_gene.txt
# OUTPUT:  carAur01.het_remove_contig.txt
# }}}
####################################################

########################################
# copy protein sequence
########################################
#******* ADD 2018-03-15 ********
cat ../maker4b/carAur01.all.maker.proteins.renamed.fasta | sed '/^>/ s/\s.*$//' > big/carAur01.maker.proteins.fasta
#cat big/carAur01.maker.proteins.fasta | sed '/^>/ s/^>/>GF|/' > big/carAur01.maker.proteins.with_sp.fasta 
cat ../maker4b/carAur01.all.maker.transcripts.renamed.fasta | sed '/^>/ s/\s.*$//' > big/carAur01.maker.transcripts.fasta
#cat carAur01.maker.transcripts.fasta | sed '/^>/ s/^>/>GF|/' > big/carAur01.maker.transcripts.with_sp.fasta 

~/my_program3/src/utility/czl_fasta_fetch -m 1 --list carAur01.gene.unmasked.tids -i big/carAur01.maker.transcripts.fasta -o big/carAur01.maker.transcripts.unmasked.fasta --cut-name -1;    mv big/carAur01.maker.transcripts.unmasked.fastafa big/carAur01.maker.transcripts.unmasked.fasta
~/my_program3/src/utility/czl_fasta_fetch -m 1 --list carAur01.gene.unmasked.tids -i big/carAur01.maker.proteins.fasta -o big/carAur01.maker.proteins.unmasked.fasta --cut-name -1;    mv big/carAur01.maker.proteins.unmasked.fastafa big/carAur01.maker.proteins.unmasked.fasta

cat carAur01.proteins.fasta.fai | awk '$2<90' > carAur01.short_protein.tid.txt
#*******************************


##########################
# get the gene set without het and not_contained
##########################
#TODO cat carAur01.gene.annot.bgp

hgFakeAgp -minContigGap=1 carAur01.withM.sm.fa carAur01.agp
hgGcPercent -wigOut -doGaps -file=stdout -win=200 -verbose=0 carAur01 carAur01.2bit | wigToBigWig stdin carAur01.chromSizes bbi/gc.bw


#######################################################
# add repeat region
#######################################################
# {{{
echo '##gff-version 3' > $asm1_dir/$asm1.gene.repeatmasker.gff
##
# old DEL 2018-03-13
# zcat $asm1_dir/../maker4/$asm1.gene.repeatmasker.gff3.gz $asm1_dir/../maker4a/$asm1.gene.repeatmasker.x1.gff.gz | grep '^[^#]'  >> $asm1_dir/$asm1.gene.repeatmasker.gff
#gzip -f $asm1_dir/$asm1.gene.repeatmasker.gff
#name=$asm1.repeatmasker
#zcat carAur01.gene.repeatmasker.gff.gz | perl -ne 'if (m/^#/) {print; next; } chomp; @t=split "\t"; if ($t[2] eq "match") {$t[2]="ncRNA";} elsif ($t[2] eq "match_part") {$t[2]="exon";} print join("\t",@t), "\n";' | gff3ToGenePred -attrsOut=$name.attr -geneNameAttr=Name stdin $name.gp
#genePredToBigGenePred $name.gp stdout | sort -k1,1 -k2,2n > $name.bgp
## remove contained/overlap region
#cat $asm1_dir/$asm1.repeatmasker.bgp | perl -e 'my $chr; $name;$b=0; $e=0; $t0; 
#while(<>) {
#    chomp;
#    my @t=split "\t",$_,-1;
#    my $chr1 = $t[0];
#    my $b1 = $t[1];
#    my $e1 = $t[2];
#    my $name1= $t[12];
#    if (defined $chr && $chr eq $chr1) {
#        my $b2 = $b1>$b ? $b1 : $b;
#        my $e2 = $e1<$e ? $e1 : $e;
#        if ($e2>$b2 && ($e2-$b2)>($e1-$b1)*0.5 && $name1 eq $name) {
#        } else {
#            if (defined $chr) { print "$t0\n"; }
#            $chr = $chr1;
#            $b   = $b1;
#            $e   = $e1;
#            $name= $name1;
#            $t0  = $_;
#        }
#    } else {
#        if (defined $chr) { print "$t0\n"; }
#        $chr = $chr1;
#        $b   = $b1;
#        $e   = $e1;
#        $name= $name1;
#        $t0  = $_;
#    }
#}
#if (defined $chr) { print "$t0\n"; }
#' > $asm1_dir/$asm1.repeatmasker.f.bgp 
#bedToBigBed -extraIndex=name,name2,geneName,geneName2 -tab -as=$kent/bigGenePred.as -type=bed12+8 $asm1_dir/$asm1.repeatmasker.f.bgp $asm1.chromSizes $asm1_dir/bbi/$asm1.repeatmasker.bb
#cat carAur01.repeatmasker.attr | awk '$2=="Name"' | cut -f 1,3 | perl -e 'my %g;
#while(<>) {
#    chomp; if (m/^#/){next;} 
#    @t = split "\t";
#    push @{$g{$t[0]}}, $t[1];
#}
#foreach my $id (sort keys(%g)) { print "$id\t", join("\t",@{$g{$id}}), "\n";}
#'> carAur01.repeatmasker.attr.Name
#cat carAur01.repeatmasker.f.bgp | cut -f 4,13 | sed -e 's/|/ /g' -e 's/species://' -e 's/genus://' > carAur01.repeatmasker.attr.Name
#ixIxx $asm1_dir/$name.attr.Name $asm1_dir/bbi/$name.ix $asm1_dir/bbi/$name.ixx

#******************************************
# add 2018-03-13
zcat ../maker4b/carAur01.all.repeatmasker.gff.gz | perl -ne 'if (m/^#/) {next; } 
chomp; @t=split "\t"; 
if ($t[2] eq "match") {
my %info;
foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa; $info{$u}=$v;}
if ($t[5]=~m/^\./) { $t[5]=0; }
if ($t[5]>1000) { $t[5]=1000; }
$info{Name}=~s/species://;
$info{Name}=~s/genus://;
$info{Name}=~s/\|/__/g;
$info{Name}=~s/%28/\(/g;
$info{Name}=~s/%29/\)/g;
$info{Name}=~s/%23/__/g;
$info{Name}=~s/%2F/__/g;
print join("\t",($t[0],$t[3]-1,$t[4],$info{Name},$t[5],$t[6])), "\n";
}' | grep -v tig00037772 | sort -k1,1 -k2,2n -k3,3nr -k4,4 | uniq > big/carAur01.repeatmasker.bed

bedToBigBed -extraIndex=name big/carAur01.repeatmasker.bed carAur01.withM.sm.fa.fai bbi/carAur01.repeatmasker.bed.bb
cat big/carAur01.repeatmasker.bed | awk -F$'\t' '{print $4"\t"$4}' > big/carAur01.repeatmasker.bed.Name
ixIxx big/carAur01.repeatmasker.bed.Name bbi/carAur01.repeatmasker.bed.bb.ix bbi/carAur01.repeatmasker.bed.bb.ixx  
#******************************************


#-----------------------------------
name=$asm1.repeatrunner
cp ../maker4a/carAur01.gene.repeatrunner.gff.gz $asm1.repeatrunner.gff.gz
zcat $name.gff.gz | perl -ne 'if (m/^#/) {print; next; } chomp; @t=split "\t"; if ($t[2] eq "protein_match") {$t[2]="ncRNA";} elsif ($t[2] eq "match_part") {$t[2]="exon";} print join("\t",@t), "\n";' | gff3ToGenePred -attrsOut=$name.attr -geneNameAttr=Name stdin $name.gp
genePredToBigGenePred $name.gp stdout | sort -k1,1 -k2,2n > $name.bgp
#-----------------------------------
zcat carAur01.gene.repeatmasker.gff.gz carAur01.gene.repeatrunner.gff.gz | grep '^[^#]' | perl -ne 'chomp; @t=split "\t"; if ($t[2] eq "match") {
     @attr=split /\s*;\s*/, $t[8];
     my $name;
     foreach my $a (@attr) {
         ($u,$v) = split /\s*=\s*/, $a;
         if ($u eq "Name") { $name=$v; }
     }
     my @out=($t[0], $t[3]-1, $t[4], $name, $t[5], $t[6]);
     print join("\t",@out), "\n";
}' | bedtools sort -i stdin | uniq > $asm1.repeat.bed
bedtools merge -i $asm1.repeat.bed -c 4 -o distinct > $asm1.repeat.bedtools_merge.bed
cat carAur01.repeat.bedtools_merge.bed | awk -v OFS=$'\t' '{gsub(/species:/,"",$4); gsub(/genus:/,"",$4); if (length($4)>255) { $4=substr($4,0,255)} print $0}' > carAur01.repeat.bedtools_merge.shorten_id.bed
bedToBigBed carAur01.repeat.bedtools_merge.shorten_id.bed carAur01.sm.fa.fai bbi/carAur01.repeat.bedtools_merge.bb
# }}}

#######################################################
# tRNA
#######################################################
cp ../miRBase_run1/goldfish.arrow.renamed.masked.tRNA.gff carAur01.tRNA.gff
perl  ~/program/biocode/gff/convert_tRNAScanSE_to_gff3.pl -i ../miRBase_run1/carAur01.tRNA.out | awk -v OFS=$'\t' '{a=$1; for (i=2;i<=8;i++) {a=a"\t"$i;} b=$9; for (i=10;i<=NF;i++) {b=b" "$i;} print a,b}' > ./big/carAur01.tRNA.gff
gff3ToGenePred -defaultCdsStatusToUnknown big/carAur01.tRNA.gff big/carAur01.tRNA.gp
genePredToBigGenePred big/carAur01.tRNA.gp stdout | sort -k1,1 -k2,2n > big/carAur01.tRNA.bgp
bedToBigBed -extraIndex=name,name2,geneName,geneName2 -tab -as=$kent/bigGenePred.as -type=bed12+8 big/carAur01.tRNA.bgp carAur01.withM.sm.fa.fai bbi/carAur01.tRNA.bb

#######################################################
# miRNA
#######################################################
f=big/carAur01.noM.miRNA
b=bbi/carAur01.noM.miRNA
cat ../miRBase_run1/ex_split/all.out.psl | grep -v tig00037772 | awk -F$'\t' -v OFS=$'\t' '{gsub(/__.*$/,"",$10); print $0}'> big/carAur01.noM.miRNA.swap.psl
cat ../miRBase_run1/ex_split/all.out.psl | grep -v tig00037772 | awk -F$'\t' -v OFS=$'\t' '{print $14"_"$16"_"$17"_"$9"_"$10}' | sed 's/__/\t/g' | sort -k1,1 > a2
#pslSwap $f.psl stdout | pslPosTarget stdin $f.swap.psl
pslRecalcMatch $f.swap.psl carAur01.2bit ~/data/miRBase/hairpin.U2T.fa  $f.swap.1.psl; mv $f.swap.1.psl $f.swap.psl;
# filter by identity  (>=97%)
pslToBigPsl -fa=$HOME/data/miRBase/hairpin.U2T.fa $f.swap.psl stdout | sort -k1,1 -k2,2n | awk -F$'\t' -v OFS=$'\t' '{print $1"_"$2"_"$3"_"$6"_"$4"\t"$0}' | sort -k1,1 > a1
join -t$'\t' -j 1 a1 a2 | cut -f 2- | awk -F$'\t' -v OFS=$'\t' '{$4=$4"__"$26; print $0;}' | cut -f 1-25 | sort -k1,1 -k2,2n >  big/carAur01.noM.miRNA.bigpsl 

cat big/carAur01.noM.miRNA.bigpsl | perl -ne 'chomp; my @t=split /\t/;
my $ne = $t[9];
if ($ne==1) { print $_, "\n"; next; }
my @sz=split /,/, $t[10];
my @bs=split /,/, $t[11];
my $m = 0;
my $b = 0;
my $e = $bs[0]+$sz[0];
my @sz1; my @bs1;
for (my $i=1; $i<$ne; $i++) {
	if ($bs[$i]-$e<30) { $e = $bs[$i]+$sz[$i]; }
	else { push @bs1, $b;  push @sz1, $e-$b;  $b=$bs[$i];  $e=$b+$sz[$i]; }
}
push @bs1, $b;  push @sz1, $e-$b;
for (my $i=0; $i<@sz1; $i++) {
	if ($sz1[$i]<30) {
		if ($i==0 && $bs1[$i+1]-$bs1[$i]>30) { $m++; }
		elsif ($i==$ne-1 && $bs1[$i]-$bs1[$i-1]>30) { $m++; }
		elsif ($bs1[$i]-$bs1[$i-1]>30 & $bs1[$i+1]-$bs1[$i]>30) { $m++; }
	}
}
if ($m==0) { print $_, "\n"; }
' > big/carAur01.noM.miRNA.f.bigpsl 

cat big/carAur01.noM.miRNA.f.bigpsl | perl -ne 'chomp; my @t=split /\t/;
my $ne = $t[9];
if ($ne==1) { print join("\t", @t[0..11]), "\n"; next; }
my @sz=split /,/, $t[10];
my @bs=split /,/, $t[11];
my $m = 0;
my $b = 0;
my $e = $bs[0]+$sz[0];
my @sz1; my @bs1;
for (my $i=1; $i<$ne; $i++) {
	if ($bs[$i]-$e<30) { $e = $bs[$i]+$sz[$i]; }
	else { push @bs1, $b;  push @sz1, $e-$b;  $b=$bs[$i];  $e=$b+$sz[$i]; }
}
push @bs1, $b;  push @sz1, $e-$b;
my @t1=@t[0..8];
$t1[9] = $#sz1+1;
$t1[10]= join(",", @sz1);
$t1[11]= join(",", @bs1);
print join("\t", @t1), "\n";
' | awk -F$'\t' -v OFS=$'\t' '{ if (chr==$1) {b1=b<$2?$2:b; e1=e<$3?e:$3; if (b1>=e1 && e1-b1<0.3*($3-$2)) {print $0; b=$2; e=$3;} } else {print $0; b=$2; e=$3; chr=$1;} }'> big/carAur01.noM.miRNA.f.bed12

bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigPsl.as $f.bigpsl carAur01.withM.sm.fa.chromSizes $b.bigpsl.bb

cat $f.swap.psl | awk -F$'\t' '($1+$3)/($2+$1+$3)>=0.95' > $f.swap.I95.psl 
pslToBigPsl -fa=$HOME/data/miRBase/hairpin.U2T.fa $f.swap.I95.psl stdout | sort -k1,1 -k2,2n | awk -F$'\t' -v OFS=$'\t' '{print $1"_"$2"_"$3"_"$6"_"$4"\t"$0}' | sort -k1,1 > a1
join -t$'\t' -j 1 a1 a2 | cut -f 2- | awk -F$'\t' -v OFS=$'\t' '{$4=$4"__"$26; print $0;}' | cut -f 1-25 | sort -k1,1 -k2,2n >  big/carAur01.noM.miRNA.I95.bigpsl 
bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as $f.I95.bigpsl carAur01.withM.sm.fa.fai $b.I95.bigpsl.bb
#cat carAur01.miRNA.bigpsl | awk -F$'\t' '$10<=2 && $14-$13>$16*0.7' | sort -k1,1 -k2,2n | bedtools cluster -i stdin | perl -e '$c=-1;$p; 
#while(<>) {
#    chomp; @t=split "\t",$_,-1; 
#    my $mr; my $s; if ($t[3]=~m/^(.*)__(.*)$/) {$mr=$1;$s=$2;} 
#    if ($s<=-20) {if ($t[25]!=$c) {if ($c!=-1) {print "$p\n";} $c=$t[25]; $p=$_;} elsif ($t[3]=~m/^dre/){$p=$_;} } 
#} print "$p\n"' > carAur01.miRNA.E2.S20.bigpsl
#cat carAur01.miRNA.E2.S20.bigpsl | awk -F$'\t' '{i=index($4,"__"); mr=substr($4,1, i-1); s=substr($4,i+2); print mr}' | sort -k1,1 | uniq -c | awk '{print $2"\t"$1}' > carAur01.miRNA.hit.hist

cat big/carAur01.noM.miRNA.I95.bigpsl | perl -ne 'chomp; my @t=split /\t/;
my $ne = $t[9];
if ($ne==1) { print $_, "\n"; next; }
my @sz=split /,/, $t[10];
my @bs=split /,/, $t[11];
my $m = 0;
my $b = 0;
my $e = $bs[0]+$sz[0];
my @sz1; my @bs1;
for (my $i=1; $i<$ne; $i++) {
	if ($bs[$i]-$e<30) { $e = $bs[$i]+$sz[$i]; }
	else { push @bs1, $b;  push @sz1, $e-$b;  $b=$bs[$i];  $e=$b+$sz[$i]; }
}
push @bs1, $b;  push @sz1, $e-$b;
for (my $i=0; $i<@sz1; $i++) {
	if ($sz1[$i]<30) {
		if ($i==0 && $bs1[$i+1]-$bs1[$i]>30) { $m++; }
		elsif ($i==$ne-1 && $bs1[$i]-$bs1[$i-1]>30) { $m++; }
		elsif ($bs1[$i]-$bs1[$i-1]>30 & $bs1[$i+1]-$bs1[$i]>30) { $m++; }
	}
}
if ($m==0) { print $_, "\n"; }
' > big/carAur01.noM.miRNA.I95.f.bigpsl 

cat big/carAur01.noM.miRNA.I95.f.bigpsl | perl -ne 'chomp; my @t=split /\t/;
my $ne = $t[9];
if ($ne==1) { print join("\t", @t[0..11]), "\n"; next; }
my @sz=split /,/, $t[10];
my @bs=split /,/, $t[11];
my $m = 0;
my $b = 0;
my $e = $bs[0]+$sz[0];
my @sz1; my @bs1;
for (my $i=1; $i<$ne; $i++) {
	if ($bs[$i]-$e<30) { $e = $bs[$i]+$sz[$i]; }
	else { push @bs1, $b;  push @sz1, $e-$b;  $b=$bs[$i];  $e=$b+$sz[$i]; }
}
push @bs1, $b;  push @sz1, $e-$b;
my @t1=@t[0..8];
$t1[9] = $#sz1+1;
$t1[10]= join(",", @sz1);
$t1[11]= join(",", @bs1);
print join("\t", @t1), "\n";
' | awk -F$'\t' -v OFS=$'\t' '{ if (chr==$1) {b1=b<$2?$2:b; e1=e<$3?e:$3; if (b1>=e1 && e1-b1<0.3*($3-$2)) {print $0; b=$2; e=$3;} } else {print $0; b=$2; e=$3; chr=$1;} }'> big/carAur01.noM.miRNA.I95.f.bed12




#######################################################
# All ncRNA
#######################################################
echo '##gff-version 3' > ./big/carAur01.ncrna.gff 
cat ../ncrna_combine/carAur01.ncrna.cpc.gff | grep -v '^#' >> ./big/carAur01.ncrna.gff 
gff3ToGenePred -attrsOut=big/carAur01.ncrna.gff.attrsOut big/carAur01.ncrna.gff big/carAur01.ncrna.gp
cat big/carAur01.ncrna.gff.attrsOut | awk -F$'\t' -v OFS=$'\t' '$2=="Name" {print $1,$3,$3}' | sort -k1,1 > big/carAur01.ncrna.gff.attrsOut.genename
genePredToBigGenePred big/carAur01.ncrna.gp stdout | sort -k4,4 | awk -F$'\t' '{print $4"\t"$0}' > a
join -t$'\t' -j 1 a big/carAur01.ncrna.gff.attrsOut.genename | awk -F$'\t' '{a=$2; for (i=3;i<=18;i++) {a=a"\t"$i;} print a"\t"$22"\t"$22"\t"$21}' > a1
join -t$'\t' -v 1 -j 1 a big/carAur01.ncrna.gff.attrsOut.genename | cut -f 2- > a2
#genePredToBigGenePred -geneNames=big/carAur01.ncrna.gff.attrsOut.genename big/carAur01.ncrna.gp stdout | sort -k1,1 -k2,2n > big/carAur01.ncrna.bgp
cat a1 a2 | sort -k1,1 -k2,2n >  big/carAur01.ncrna.bgp
bedToBigBed -extraIndex=name,name2,geneName,geneName2 -tab -as=$kent/bigGenePred.as -type=bed12+8 big/carAur01.ncrna.bgp $asm1.chromSizes bbi/carAur01.ncrna.bgp.bb
cat big/carAur01.ncrna.gff | perl -ne '
if (m/^#/ || m/^\s*$/) { next; } chomp; my @t=split "\t",$_,-1;
my %info; foreach my $aa (split /\s*;\s*/,$t[8]) {my ($u,$v)=split /\s*=\s*/,$aa,2; $info{$u}=$v;}
if ($t[2]=~m/exon/) { print join("\t", ($t[0], $t[3]-1, $t[4], "$info{Parent}:$info{exon_number}", $t[5], $t[6]) ), "\n" }' | sort -k1,1 -k2,2n > big/carAur01.ncrna.exon.bed

#######################################################
# Rfam
#######################################################
f=big/carAur01.noM.rfam
b=bbi/carAur01.noM.rfam
cat ../miRBase_run1/carAur01.rfam.cmscan.f.bgp | grep -v tig00037772 > $f.bgp
bedToBigBed -extraIndex=name,name2 -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as $f.bgp carAur01.withM.sm.fa.fai $b.bb
#cat $fn.gff.attrs | awk '$2=="ID"' | cut -f 1,3 > $fn.attrs.ID
#ixIxx $fn.attrs.ID bbi/$asm1.rfam.ix bbi/$asm1.rfam.ixx

################################################
#  add Interpro and GO
################################################
#****************************************
# convert
#****************************************
cat ../maker4b/carAur01.ips.gff | awk -F$'\t' -v OFS=$'\t' '$1~/^#/ {print} $3=="match" {$3="mRNA";print} $3=="match_part" {$3="exon"; print}' | gff3ToGenePred -useName -defaultCdsStatusToUnknown stdin stdout | grep -v tig00037772 > big/carAur01.noM.ips.gp
genePredToBigGenePred big/carAur01.noM.ips.gp stdout  | sort -k1,1 -k2,2n > big/carAur01.noM.ips.bgp
bedToBigBed -extraIndex=name,name2,geneName,geneName2 -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as big/carAur01.noM.ips.bgp carAur01.withM.sm.fa.fai big/carAur01.noM.ips.bgp.bb
#
cat ../maker4b/carAur01.ips.gff | awk -F$'\t' -v OFS=$'\t' '$1~/^#/ {print} $3=="match" {$3="mRNA";print} $3=="match_part" {$3="exon"; print}' | gff3ToGenePred -rnaNameAttr=Ontology_term -defaultCdsStatusToUnknown stdin stdout | grep -v tig00037772 > big/carAur01.noM.GO.gp
genePredToBigGenePred big/carAur01.noM.GO.gp stdout  | sort -k1,1 -k2,2n > big/carAur01.noM.GO.bgp
cat big/carAur01.noM.GO.bgp | awk -F$'\t' '{print $18"\t"$0}' | sort -k1,1 > a;
join -t$'\t' -j1 a ~/data/GO/go.tsv | cut -f 2-21,23,24 | awk -F$'\t' -v OFS=$'\t' '{gsub(/ /,"_",$22); $18=$18":"$21":"$22; print $0}' | cut -f 1-20  | sort -k1,1 -k2,2n > big/carAur01.noM.GO.rename.bgp 

bedToBigBed -extraIndex=name,geneName -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as big/carAur01.noM.GO.rename.bgp carAur01.withM.sm.fa.fai big/carAur01.noM.GO.bgp.bb
#****************************************

#****************************************
# add 2018-04-04
# Interpro + GO + Pathway
#****************************************
#{{{
mkdir big/ips
/home/chenz11/my_program3/src/utility/czl_GO_to_tab.pl -i ~/data/annot_db2/go/go.obo -o GO. ;
mv GO.molecular_function.txt GO_MF.txt
mv GO.biological_process.txt GO_BP.txt
mv GO.cellular_component.txt GO_CC.txt
for d in MF BP CC
do
    for i in 1 2 3
    do
        cat GO_$d.txt | awk -F$'\t' -v OFS=$'\t' '$3=='$i' {print $1,$2}' | sort -k1,1 | uniq > GO_$d.lv$i.txt;
    done
done
#
cat ../maker4b/carAur01.ips.rename.tsv | sort -k1,1 > ./big/ips/carAur01.ips.tsv
join -t$'\t' -j 1 -v 1 big/ips/carAur01.ips.tsv carAur01.gene.masked.tids > big/ips/carAur01.ips.unmasked.tsv 
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i big/ips/carAur01.ips.unmasked.tsv -o big/ips/carAur01.ips.unmasked. -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt --max-go-lv 6
cat big/carAur01.gene.gff.gtpnnbbcc | awk '$3!="" && $3!="." {print $3"\t"$1}' | sort -k1,1 > a;
join -t$'\t' -j 1 a  big/ips/carAur01.ips.unmasked.tsv | cut -f 2- | sort -k1,1 > big/ips/carAur01.gene.ips.unmasked.tsv 
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i big/ips/carAur01.gene.ips.unmasked.tsv -o big/ips/carAur01.gene.ips.unmasked. -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt
cat big/ips/carAur01.gene.ips.unmasked.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > big/ips/carAur01.gene.ips.unmasked.id_to_Interpro.txt
cat big/ips/carAur01.gene.ips.unmasked.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > big/ips/carAur01.gene.ips.unmasked.id_to_Pfam.txt
cat big/ips/carAur01.gene.ips.unmasked.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > big/ips/carAur01.gene.ips.unmasked.id_to_Gene3D.txt
mkdir big/ips/compare
mkdir big/ips/compare/merged
for f in Reactome KEGG Interpro Pfam GO2_CC.lv2 GO2_BP.lv2 GO2_MF.lv2 GO2_CC.lv3 GO2_BP.lv3 GO2_MF.lv3 GO2_CC.lv4 GO2_BP.lv4 GO2_MF.lv4 GO2_CC.lv5 GO2_BP.lv5 GO2_MF.lv5 GO2_CC.lv6 GO2_BP.lv6 GO2_MF.lv6
do
	echo $f
	for sp1 in ZF GC CC GF
	do
	for sp2 in ZF GC CC GF
	do
		if [ "$sp1" == "$sp2" ]; then continue; fi
		if ! [ -d "big/ips/compare/$sp1.$sp2" ]; then mkdir "big/ips/compare/$sp1.$sp2"; fi
		Rscript ~/my_program3/src/annot_genome/czl_cmp_two_set.Rscript big/ips/$sp1/gene.ips.id_to_$f.txt big/ips/$sp2/gene.ips.id_to_$f.txt big/ips/compare/$sp1.$sp2/$f.enrich.txt
	done
	done
    cat big/ips/compare/ZF.GC/$f.enrich.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' -v OFS=$'\t' '{if ($7==0) $7=0.001; else if ($7~/^Inf/) $7=1000; print $0}'> a1
    cat big/ips/compare/ZF.GF/$f.enrich.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' -v OFS=$'\t' '{if ($7==0) $7=0.001; else if ($7~/^Inf/) $7=1000; print $0}' > a2
    cat big/ips/compare/ZF.CC/$f.enrich.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' -v OFS=$'\t' '{if ($7==0) $7=0.001; else if ($7~/^Inf/) $7=1000; print $0}' > a3
    cat big/ips/compare/GC.GF/$f.enrich.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' -v OFS=$'\t' '{if ($7==0) $7=0.001; else if ($7~/^Inf/) $7=1000; print $0}' > a4
    cat big/ips/compare/GC.CC/$f.enrich.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' -v OFS=$'\t' '{if ($7==0) $7=0.001; else if ($7~/^Inf/) $7=1000; print $0}' > a5
    cat big/ips/compare/GF.CC/$f.enrich.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' -v OFS=$'\t' '{if ($7==0) $7=0.001; else if ($7~/^Inf/) $7=1000; print $0}' > a6
    perl -e '$sp1="ZF"; $sp2="GC"; print join("\t", ("ID", "name", "$sp1.n1", "$sp2.n1", "$sp1.n2", "$sp2.n2", "${sp1}_${sp2}.odd", "${sp1}_${sp2}.p", "${sp1}_${sp2}.FDR", "$sp1.perc1", "$sp2.perc1"));'> h1
    perl -e '$sp1="ZF"; $sp2="GF"; print join("\t", ("ID", "name", "$sp1.n1", "$sp2.n1", "$sp1.n2", "$sp2.n2", "${sp1}_${sp2}.odd", "${sp1}_${sp2}.p", "${sp1}_${sp2}.FDR", "$sp1.perc1", "$sp2.perc1"));'> h2
    perl -e '$sp1="ZF"; $sp2="CC"; print join("\t", ("ID", "name", "$sp1.n1", "$sp2.n1", "$sp1.n2", "$sp2.n2", "${sp1}_${sp2}.odd", "${sp1}_${sp2}.p", "${sp1}_${sp2}.FDR", "$sp1.perc1", "$sp2.perc1"));'> h3
    perl -e '$sp1="GC"; $sp2="GF"; print join("\t", ("ID", "name", "$sp1.n1", "$sp2.n1", "$sp1.n2", "$sp2.n2", "${sp1}_${sp2}.odd", "${sp1}_${sp2}.p", "${sp1}_${sp2}.FDR", "$sp1.perc1", "$sp2.perc1"));'> h4
    perl -e '$sp1="GC"; $sp2="CC"; print join("\t", ("ID", "name", "$sp1.n1", "$sp2.n1", "$sp1.n2", "$sp2.n2", "${sp1}_${sp2}.odd", "${sp1}_${sp2}.p", "${sp1}_${sp2}.FDR", "$sp1.perc1", "$sp2.perc1"));'> h5
    perl -e '$sp1="GF"; $sp2="CC"; print join("\t", ("ID", "name", "$sp1.n1", "$sp2.n1", "$sp1.n2", "$sp2.n2", "${sp1}_${sp2}.odd", "${sp1}_${sp2}.p", "${sp1}_${sp2}.FDR", "$sp1.perc1", "$sp2.perc1"));'> h6
	cp a1 b;
	cp h1 h;
	for i in `seq 2 6`
	do
		join -t $'\t' -j 1 b a$i > b1; mv b1 b;
		join -t $'\t' -j 1 h h$i > h1; mv h1 h;
	done
	cat h b > big/ips/compare/merged/$f.enrich.txt
#    head -n 1 big/ips/compare/GF.CC.$f.enrich.txt | sort -k1,1 > h1
#    head -n 1 big/ips/compare/GF.ZF.$f.enrich.txt | sort -k1,1 > h2
#    join -t $'\t' -j 1 h1 h2 > big/ips/compare/GF.CC.ZF.$f.enrich.merged.txt 
#    cat big/ips/compare/GF.CC.$f.enrich.txt | tail -n +2 | sort -k1,1 > a1
#    cat big/ips/compare/GF.ZF.$f.enrich.txt | tail -n +2 | sort -k1,1 > a2
#    join -t $'\t' -j 1 a1 a2 >> big/ips/compare/GF.CC.ZF.$f.enrich.merged.txt 
done
#
for f in CC BP MF
do
    for i in 1 2 3
    do 
        n=`cat big/ips/carAur01.gene.ips.unmasked.id_to_GO2_${f}.lv$i.txt | cut -f 1 | sort | uniq | wc -l`
        cat big/ips/carAur01.gene.ips.unmasked.id_to_GO2_${f}.lv$i.txt | awk -F$'\t' -v OFS=$'\t' '{print $2"___"$3}' | sed 's/ /__/g' | sort -k1,1 | uniq -c | awk '{print $2"\t"$1}' | awk -F$'\t' -v OFS=$'\t' -v n=$n '{print $1,$2,n,$3,n-$3,int($3*1000000/n)/10000}' > big/ips/carAur01.gene.ips.unmasked.GO_$f.lv$i.counts.txt
    done
done

# merge all levels of GO
for gs in GO2_CC GO2_BP GO2_MF
do
cat big/ips/compare/merged/$gs.lv$i.enrich.txt | head -n 1 | awk '{print $0"\tGO_level"}' > big/ips/compare/merged/$gs.enrich.txt
for i in `seq 2 6`
do
	cat big/ips/compare/merged/$gs.lv$i.enrich.txt | tail -n +2 | awk '{print $0"\t"'$i'}' >> big/ips/compare/merged/$gs.enrich.txt
done
done
for gs in GO2_CC GO2_BP GO2_MF KEGG Reactome Interpro Pfam
do
	cat big/ips/compare/merged/$gs.enrich.txt | awk -F$'\t' 'NR==1 {print} NR>1 { n=0; for (i=0;i<6;i++) {j1=8+i*10; j2=9+i*10; if ($j1<0.01 && $j2<0.1) {n++;} } if (n>0) {print}}' > big/ips/compare/merged/$gs.enrich.f.txt
done
#}}}
#****************************************

mkdir big/ips_run2
mkdir big/ips_run2/protein.ips.parse
cat ../maker4b/interProScan5_run2/carAur01.p.ips/carAur01.all.maker.proteins.renamed.fasta.*.tsv | sort -k1,1 > big/ips_run2/protein.ips.parse/protein.ips.tsv
cd big/ips_run2/
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i protein.ips.tsv -o protein.ips.parse/protein.ips. -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list
cat protein.ips.tsv | sed 's/_R[0-9]\+\t/\t/' > gene.ips.tsv
perl ~/my_program3/src/annot_genome/czl_parse_ips_tsv.pl -i gene.ips.tsv -o protein.ips.parse/gene.ips. -go ~/data/annot_db2/go/go.obo -kegg ~/data/annot_db2/kegg/map/pathway.list -reactome ~/data/reactome/ReactomePathways.txt
cat gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$12!="" {print $1,$12,$13}' | sort -k1,1 -k2,2 | uniq > protein.ips.parse/gene.ips.id_to_Interpro.txt
cat gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Pfam" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > protein.ips.parse/gene.ips.id_to_Pfam.txt
cat gene.ips.tsv | awk -F$'\t' -v OFS=$'\t' '$4=="Gene3D" {print $1,$5,$6}' | sort -k1,1 -k2,2 | uniq > protein.ips.parse/gene.ips.id_to_Gene3D.txt




#######################################################
# process Trinity-GG exonerate to psl
#######################################################
#{{{
for sm in `cat ../Trinity_GG_run01/sample`
do
f=big/est2genome/$sm.trinity_GG
b=bbi/est2genome/$sm.trinity_GG
cat ../Trinity_GG_run01/$sm.trinity/Trinity-GG.carAur01.ex2/all.out.f.match.gff3 | grep -v tig00037772 | gff3ToPsl carAur01.withM.sm.fa.fai ../Trinity_GG_run01/$sm.trinity/Trinity-GG.short_name.fasta.fai stdin big/est2genome/$sm.trinity_GG.psl
pslSwap big/est2genome/$sm.trinity_GG.psl stdout | pslPosTarget stdin stdout | awk -F$'\t' -v OFS=$'\t' '{gsub(/\.[0-9]+$/,"",$10); print $0}' > big/est2genome/$sm.trinity_GG.swap.psl 
pslRecalcMatch $f.swap.psl carAur01.2bit ../Trinity_GG_run01/$sm.trinity/Trinity-GG.short_name.fasta  $f.swap.1.psl; mv $f.swap.1.psl $f.swap.psl;
# filter by identity  (>=97%)
cat $f.swap.psl | awk -F$'\t' '($1+$3)/($2+$1+$3)>=0.97' > $f.swap.I97.psl 
pslToBigPsl -fa=../Trinity_GG_run01/$sm.trinity/Trinity-GG.short_name.fasta $f.swap.I97.psl stdout | sort -k1,1 -k2,2n > $f.bigpsl
bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as $f.bigpsl carAur01.withM.sm.fa.fai $b.bigpsl.bb
done
#}}}
#######################################################

#######################################################
# add RNA-seq EST Alignment
#######################################################
#***************************
#  delete 2018-03-17
#***************************
#echo '##gff-version 3' > $asm1_dir/$asm1.gene.est2genome.gff
#zcat $asm1_dir/../maker4/$asm1.gene.est2genome.gff3.gz $asm1_dir/../maker4a/$asm1.gene.est2genome.x1.gff.gz | grep '^[^#]' >> $asm1_dir/$asm1.gene.est2genome.gff
#gzip -f $asm1_dir/$asm1.gene.est2genome.gff
#zcat carAur01.gene.est2genome.gff.gz | perl -ne 'if (m/^#/) {print; next; } chomp; @t=split "\t"; if ($t[2] eq "expressed_sequence_match") {$t[2]="match";} print join("\t",@t), "\n";' | sed 's/Length=/length=/' | perl ~/my_program3/src/utility/czl_gff_match_part_to_match.pl -i - -o $asm1_dir/$asm1.gene.est2genome.1.gff.gz
#zcat $asm1_dir/$asm1.gene.est2genome.1.gff.gz | gff3ToPsl "$asm1_dir/$asm1.chromSizes" "/home/chenz11/data/goldfish/RNA_seq/trinity_run1/merge9.cdhit.short_name.fa.fai" stdin stdout | gzip -c > $asm1.est2genome.psl.gz
#zcat $asm1.est2genome.psl.gz | pslSwap stdin stdout | pslPosTarget stdin $asm1.est2genome.psl.tmp
#mv $asm1.est2genome.psl.tmp $asm1.est2genome.psl; gzip -f $asm1.est2genome.psl
#zcat $asm1_dir/$asm1.est2genome.psl.gz | pslToBigPsl stdin stdout | sort -k1,1 -k2,2n | gzip -c > $asm1.est2genome.bigpsl.gz
#zcat $asm1_dir/$asm1.est2genome.bigpsl.gz > $asm1_dir/$asm1.est2genome.bigpsl
#bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=$kent/bigPsl.as "$asm1_dir/$asm1.est2genome.bigpsl" "$asm1_dir/$asm1.chromSizes" "$asm1_dir/bbi/$asm1.est2genome.bb"
#***************************

#***************************
#  add 2018-03-17
#***************************
cat ../maker4b/carAur01.all.est2genome.gff | grep -v tig00037772 | sed 's/expressed_sequence_match/match/' | sed 's/Length=/length=/g' | perl ~/my_program3/src/utility/czl_gff_match_part_to_match.pl -i - -o big/carAur01.est2genome.gff
gff3ToPsl carAur01.withM.sm.fa.fai ~/data/goldfish/RNA_seq/trinity_run1/merge9.cdhit.fai big/carAur01.est2genome.gff  big/carAur01.est2genome.psl
pslSwap big/carAur01.est2genome.psl stdout | pslPosTarget stdin stdout | awk -F$'\t' -v OFS=$'\t' '{gsub(/\.[0-9]+$/,"",$10); print $0}' > big/carAur01.est2genome.swap.psl;
pslRecalcMatch big/carAur01.est2genome.swap.psl carAur01.2bit  $HOME/data/goldfish/RNA_seq/trinity_run1/merge9.cdhit  big/carAur01.est2genome.swap.1.psl; mv big/carAur01.est2genome.swap.1.psl big/carAur01.est2genome.swap.psl;
# filter by identity  (>=97%)
cat big/carAur01.est2genome.swap.psl | awk -F$'\t' '($1+$3)/($2+$1+$3)>=0.97' > big/carAur01.est2genome.swap.I97.psl 
pslToBigPsl -fa=$HOME/data/goldfish/RNA_seq/trinity_run1/merge9.cdhit big/carAur01.est2genome.swap.I97.psl stdout | sort -k1,1 -k2,2n > big/carAur01.est2genome.bigpsl
cat big/carAur01.est2genome.bigpsl | sed -e 's/0[136]_TRI/Brain_TRI/g' -e 's/02_TRI/Eye_TRI/g' -e 's/0[47]_TRI/Gill_TRI/g' -e 's/05_TRI/Bone_TRI/g' -e 's/08_TRI/Heart_TRI/g' -e 's/09_TRI/Muscle_TRI/g' -e 's/10_TRI/TailFin_TRI/g' > big/carAur01.est2genome.rename.bigpsl 
bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as big/carAur01.est2genome.rename.bigpsl carAur01.withM.sm.fa.fai bbi/carAur01.est2genome.bigpsl.bb
#***************************

#######################################################
# add cdna Alignment
#######################################################
#********************
# delete 2018-06-15
#********************
#{{{
#echo '##gff-version 3' > $asm1_dir/$asm1.gene.cdna2genome.gff
#zcat $asm1_dir/../maker4/$asm1.gene.cdna2genome.gff3.gz $asm1_dir/../maker4a/$asm1.gene.cdna2genome.x1.gff.gz | grep '^[^#]' >> $asm1_dir/$asm1.gene.cdna2genome.gff
#gzip -f $asm1_dir/$asm1.gene.cdna2genome.gff
#zcat carAur01.gene.cdna2genome.gff.gz | perl -ne 'if (m/^#/) {print; next; } chomp; @t=split "\t"; if ($t[2] eq "expressed_sequence_match") {$t[2]="match";} print join("\t",@t), "\n";' | sed 's/Length=/length=/' | perl ~/my_program3/src/utility/czl_gff_match_part_to_match.pl -i - -o $asm1_dir/$asm1.gene.cdna2genome.1.gff.gz
#zcat $asm1_dir/$asm1.gene.cdna2genome.1.gff.gz | gff3ToPsl "$asm1_dir/$asm1.chromSizes" "$HOME/data/goldfish/other/for_gene_annotation/rna.short_name.fasta.fai" stdin stdout | gzip -c > $asm1.cdna2genome.psl.gz
#zcat $asm1.cdna2genome.psl.gz | pslSwap stdin stdout | pslPosTarget stdin $asm1.cdna2genome.psl.tmp
#mv $asm1.cdna2genome.psl.tmp $asm1.cdna2genome.psl; gzip -f $asm1.cdna2genome.psl
#zcat $asm1_dir/$asm1.cdna2genome.psl.gz | pslToBigPsl stdin stdout | sort -k1,1 -k2,2n > $asm1.cdna2genome.bigpsl
#bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=$kent/bigPsl.as "$asm1_dir/$asm1.cdna2genome.bigpsl" "$asm1_dir/$asm1.chromSizes" "$asm1_dir/bbi/$asm1.cdna2genome.bb"
#cp $asm1_dir/bbi/$asm1.cdna2genome.bb ~/data/datashare/fishCarAur1/$asm1/bbi/
#}}}
#*******************

#***************************
#  add 2018-03-16
#***************************
cat ../maker4b/carAur01.all.cdna2genome.ZF.gff | grep -v tig00037772 | sed 's/expressed_sequence_match/match/' | sed 's/Length=/length=/g' | perl ~/my_program3/src/utility/czl_gff_match_part_to_match.pl -i - -o big/carAur01.cdna2genome.ZF.gff
gff3ToPsl carAur01.withM.sm.fa.fai ~/data/goldfish/other/for_gene_annotation/rna.short_name.fasta.chromSizes.2  big/carAur01.cdna2genome.ZF.gff  big/carAur01.cdna2genome.ZF.psl
pslSwap big/carAur01.cdna2genome.ZF.psl stdout | pslPosTarget stdin stdout | awk -F$'\t' -v OFS=$'\t' '{gsub(/\.[0-9]+$/,"",$10); print $0}' > big/carAur01.cdna2genome.ZF.swap.psl;
pslToBigPsl -cds=$HOME/data/ensembl85/cds_for_bigpsl/Danio_rerio.GRCz10.85.cds -fa=/home/chenz11/data/ensembl85/ens85.rna.short_name.no_version.fa big/carAur01.cdna2genome.ZF.swap.psl stdout | sort -k1,1 -k2,2n > big/carAur01.cdna2genome.ZF.bigpsl
bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as big/carAur01.cdna2genome.ZF.bigpsl carAur01.withM.sm.fa.fai bbi/carAur01.cdna2genome.ZF.bigpsl.bb
#
f=big/carAur01.cdna2genome
b=bbi/carAur01.cdna2genome
cat ../maker4b/carAur01.all.cdna2genome.gff | grep -v tig00037772 | sed 's/expressed_sequence_match/match/' | sed 's/Length=/length=/g' | perl ~/my_program3/src/utility/czl_gff_match_part_to_match.pl -i - -o big/carAur01.cdna2genome.gff
gff3ToPsl carAur01.withM.sm.fa.fai /home/chenz11/data/goldfish/other/for_gene_annotation/rna.short_name.fasta.chromSizes.2  big/carAur01.cdna2genome.gff  big/carAur01.cdna2genome.psl
#pslSwap big/carAur01.cdna2genome.psl stdout | pslPosTarget stdin stdout | awk -F$'\t' -v OFS=$'\t' '{gsub(/\.[0-9]+$/,"",$10); print $0}' > big/carAur01.cdna2genome.swap.psl;
pslSwap big/carAur01.cdna2genome.psl stdout | pslPosTarget stdin stdout > big/carAur01.cdna2genome.swap.psl;
pslRecalcMatch $f.swap.psl carAur01.2bit /home/chenz11/data/goldfish/other/for_gene_annotation/rna.short_name.fasta $f.swap.1.psl; mv $f.swap.1.psl $f.swap.psl;
pslToBigPsl -fa=$HOME/data/goldfish/other/for_gene_annotation/rna.short_name.fasta big/carAur01.cdna2genome.swap.psl stdout | sort -k1,1 -k2,2n > big/carAur01.cdna2genome.bigpsl
bedToBigBed -type=bed12+13 -extraIndex=name -tab -as=$kent/bigPsl.as big/carAur01.cdna2genome.ZF.bigpsl carAur01.withM.sm.fa.fai bbi/carAur01.cdna2genome.ZF.bigpsl.bb
#***************************


#######################################################
# add phastcon peaks and region
#######################################################

#*************************************
# add 2018-03-26
#*************************************
#{{{
# net.chain phastcons from:
#   chain_net1/multi_run2/ZF_GC_CC_GF.netSubset.out/all.orig.roast_multic.GF.phastcons.bb => ZF_GC_CC_GF.net.GF.bb
#   chain_net1/multi_run2/GC_CC_GF.netSubset.out/all.orig.roast_multic.GF.phastcons.bb => GC_CC_GF.net.GF.bb
#   chain_net1/multi_run2/ZF_GC_GF.netSubset.out/all.orig.roast_multic.GF.phastcons.bb => ZF_GC_GF.net.GF.bb
#}}}

#************************************
# del 2018-03-26
#************************************
# {{{
#cp $asm1_dir/../chain_net1/multiz4/gcarpF.carpNG.carAur/roast.mz.carAur.phastcons.bed $asm1_dir/$asm1.GC_GF_CC.roast.mz.phastcons.bed
#bedToBigBed $asm1_dir/$asm1.ZF_GC_GF_CC.roast.mz.phastcons.bed $asm1_dir/$asm1.chromSizes $asm1_dir/bbi/$asm1.GC_GF_CC.roast.mz.phastcons.bb
#cp $asm1_dir/../chain_net1/multiz4/gcarpF.carpNG.carAur/tba.mz.carAur.phastcons.bed $asm1_dir/$asm1.GC_GF_CC.tba.mz.phastcons.bed
#bedToBigBed $asm1_dir/$asm1.ZF_GC_GF_CC.tba.mz.phastcons.bed $asm1_dir/$asm1.chromSizes $asm1_dir/bbi/$asm1.GC_GF_CC.tba.mz.phastcons.bb
#
#cat $asm1_dir/../chain_net1/multiz4/danRer10.gcarpF.carAur.carpNG/roast.mz.carAur.phastcons.bed | sort -k1,1 -k2,2n > $asm1_dir/$asm1.ZF_GC_GF_CC.roast.mz.phastcons.bed
#bedToBigBed $asm1_dir/$asm1.ZF_GC_GF_CC.roast.mz.phastcons.bed $asm1_dir/$asm1.chromSizes $asm1_dir/bbi/$asm1.ZF_GC_GF_CC.roast.mz.phastcons.bb
#cat $asm1_dir/../chain_net1/multiz4/danRer10.gcarpF.carAur.carpNG/tba.mz.carAur.phastcons.bed | sort -k1,1 -k2,2n > $asm1_dir/$asm1.ZF_GC_GF_CC.tba.mz.phastcons.bed
#bedToBigBed $asm1_dir/$asm1.ZF_GC_GF_CC.tba.mz.phastcons.bed $asm1_dir/$asm1.chromSizes $asm1_dir/bbi/$asm1.ZF_GC_GF_CC.tba.mz.phastcons.bb
#
#cp $asm1_dir/../chain_net1/multiz4/danRer10.gcarpF.carAur.carpNG/roast.mz.carAur.phastcons.bw $asm1_dir/bbi/$asm1.ZF_GC_GF_CC.roast.mz.phastcons.bw
#cp $asm1_dir/../chain_net1/multiz4/danRer10.gcarpF.carAur.carpNG/tba.mz.carAur.phastcons.bw $asm1_dir/bbi/$asm1.ZF_GC_GF_CC.tba.mz.phastcons.bw
#cp $asm1_dir/bbi/$asm1.ZF_GC_GF_CC.tba.mz.phastcons.bw ~/data/datashare/fishCarAur1/$asm1/bbi/
#cp $asm1_dir/bbi/$asm1.ZF_GC_GF_CC.roast.mz.phastcons.bw ~/data/datashare/fishCarAur1/$asm1/bbi/
#}}}

#######################################################
# chain
#######################################################
~/my_program3/src/utility/chain_to_bb.sh $asm1_dir/../chain_net1/carAur.vs.danRer10/all.target.net.chain.gz $asm1_dir/$asm1.chromSizes  $asm1_dir/bbi/vs.danRer10.net $kent
~/my_program3/src/utility/chain_to_bb.sh $asm1_dir/../chain_net1/carAur.vs.danRer10/all.prenet.chain.gz $asm1_dir/$asm1.chromSizes  $asm1_dir/bbi/vs.danRer10.prenet $kent
~/my_program3/src/utility/chain_to_bb.sh $asm1_dir/../chain_net1/carAur.vs.carAur/all.target.net.fself.chain $asm1_dir/$asm1.chromSizes  $asm1_dir/bbi/vs.carAur01.net $kent
~/my_program3/src/utility/chain_to_bb.sh $asm1_dir/../chain_net1/carAur.vs.carAur/all.prenet.fself.chain.gz $asm1_dir/$asm1.chromSizes  $asm1_dir/bbi/vs.carAur01.prenet $kent


#######################################################
# stats
#######################################################
>$asm1.stats
cat $asm1.chromSizes | awk -v l=0 '{l=l+$2} END{print "Total Genome Size (bp) = "l}' >> $asm1.stats
sz=`cat $asm1.stats | grep "Total Genome Size" | head -n 1 | sed 's/[^0-9]//g'`
cat $asm1.repeat.bedtools_merge.bed | awk -v l=0 -v sz=$sz '{l=l+($3-$2)} END{print "Total Repeat (bp) = "l" ("l*100/sz"%)"}' >> $asm1.stats


###########################
# fetch transcript ids (coding)
###########################
cat carAur01.transcripts.fasta |  grep '^>' | sed -e 's/^>//' -e 's/\s.*$//' | sort | uniq > carAur01.transcripts.fasta.ids

###########################
# GO stat
###########################
# gene transcript protein  go_domain goterm  goname
cat ../maker4a/carAur01.interProScan5.add_name.tsv | gzip -c > carAur01.interProScan5.add_name.tsv.gz
zcat carAur01.interProScan5.add_name.tsv.gz | awk -F$'\t' '$14~/GO/ {print $1"\t"$14}' | perl -ne '
chomp; @t=split "\t"; @go=split /\|/, $t[1]; $g=$t[0]; $g=~s/-mRNA.*$//;
foreach $go (@go) {@a=split ":",$go; print "$g\t$t[0]\t$t[0]\t$a[2]\t$a[0]:$a[1]\t$a[3]\n"}' |uniq > carAur01.interProScan5.GO.txt
cat carAur01.interProScan5.GO.txt | cut -f 1 | sort | uniq | wc -l | cut -d" " -f 1 > carAur01.interProScan5.GO.Ngene
cat carAur01.interProScan5.GO.txt | sort -t$'\t' -k4,4 -k5,5 | awk -F$'\t' '{print $4"__"$5"__"$6}' | uniq -c | awk '{print $2"\t"$1}' | sed 's/__/\t/g' | awk -F$'\t' -v n=`cat carAur01.interProScan5.GO.Ngene` '{print $0"\t"n-$NF}'> carAur01.interProScan5.GO.count_by_GO.txt
###########################

###########################
# Interpro stat
###########################
# gene transcript protein  go_domain goterm  goname
zcat carAur01.interProScan5.add_name.tsv.gz | awk -F$'\t' '$12~/.../ {print $1"\t"$12"\t"$13}' | uniq | perl -ne '
chomp; @t=split "\t"; $g=$t[0]; $g=~s/-mRNA.*$//;
print "$g\t$t[0]\t$t[0]\t$t[1]\t$t[2]\n";' | sort -t$'\t' -k1,1 -k2,2 | uniq > carAur01.interProScan5.interPro.txt
cat carAur01.interProScan5.interPro.txt | cut -f 1 | sort | uniq | wc -l | cut -d" " -f 1 > carAur01.interProScan5.interPro.Ngene
cat carAur01.interProScan5.interPro.txt | sort -t$'\t' -k4,4 -k5,5 | awk -F$'\t' '{print $4"__"$5}' | uniq -c | awk '{print $2"\t"$1}' | sed 's/__/\t/g' | awk -F$'\t' -v n=`cat carAur01.interProScan5.interPro.Ngene` '{print $0"\t"n-$NF}'> carAur01.interProScan5.interPro.count_by_interPro.txt
#########################


########### copy to UCSC #############
files=""
files=$files  carAur01.2bit  carAur01.agp
files=$files  bbi/carAur01.repeatmasker.bed.bb  bbi/carAur01.repeatmasker.bed.bb.ix  bbi/carAur01.repeatmasker.bed.bb.ixx
files=$files  bbi/gc.bw
files=$files  bbi/carAur01.repeat.bedtools_merge.bb
files=$files  bbi/carAur01.contigs.readToTig.bb  bbi/carAur01.contigs.readToTig.bw
files=$files  bbi/carAur01.noM.contigs.readToTig.bb  bbi/carAur01.noM.contigs.readToTig.bw
files=$files  bbi/carAur01.noM.gene.bb
files=$files  bbi/carAur01.cdna2genome.ZF.bigpsl.bb
files=$files  bbi/phastcons/ZF_GC_CC_GF.net.GF.bb  bbi/phastcons/ZF_GC_CC_GF.net.GF.bw
files=$files  bbi/phastcons/GC_CC_GF.net.GF.bb  bbi/phastcons/GC_CC_GF.net.GF.bw
files=$files  bbi/phastcons/ZF_GC_GF.net.GF.bb  bbi/phastcons/ZF_GC_GF.net.GF.bw
for f in $files
do
    cp $f ~/data/datashare/fishCarAur1/carAur01/bbi/
done
