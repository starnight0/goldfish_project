datadir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/maker4a
cwd=`pwd`
# cat two contigs to maker4 results
# tig00002911_arrow
# tig00214208_arrow
cd $datadir
ctg=tig00002911_arrow
f=`cat contig2.success | grep $ctg | cut -f 2`
types=`ls goldfish.arrow.renamed.maker.output/$f/$ctg.*.fasta | sed -e 's/^.*'$ctg'\/\+//'`
echo $types
for t1 in $types
do
  if ! [ -f carAur01.all.$t1 ] 
  then
    cp ../maker4/carAur01.all.$t1 ./
    for ctg in tig00002911_arrow tig00214208_arrow
    do
      f=`cat contig2.success | grep $ctg | cut -f 2`
      cat goldfish.arrow.renamed.maker.output/$f/$ctg.$t1 > carAur01.all.$t1
    done
  fi
done
if ! [ -f "carAur01.gene.gff.gz" ] && ! [ -f  "carAur01.gene.gff" ]
then
  echo aaaa
  exit 0
  if [ -f ../maker4/carAur01.gene.gff3.gz ]; then zcat ../maker4/carAur01.gene.gff3.gz | awk 'BEGIN {r=0} { if (r==0 && $0~/^##.*FASTA/) {r=1;} if (r==0) {print} }' >> carAur01.gene.gff; fi
  for ctg in tig00002911_arrow tig00214208_arrow
  do
    t="goldfish.arrow.renamed.maker.output/$f/$ctg.gff"
    cat $t | awk 'BEGIN {r=0} { if (r==0 && $0~/^##.*FASTA/) {r=1;} if (r==0) {print} }' >> carAur01.gene.gff
  done
fi

echo '##gff-version 3' > carAur01.gene.x1.gff
# don't run this again because ./maker4/goldfish.arrow.renamed.maker.output/goldfish.arrow.renamed_datastore/FB/57/tig00004487_arrow/tig00004487_arrow.gff is removed
# cat   goldfish.arrow.renamed.maker.output/goldfish.arrow.renamed_datastore/B3/74/tig00002911_arrow/tig00002911_arrow.gff   goldfish.arrow.renamed.maker.output/goldfish.arrow.renamed_datastore/81/5B/tig00214208_arrow/tig00214208_arrow.gff ../maker4/goldfish.arrow.renamed.maker.output/goldfish.arrow.renamed_datastore/FB/57/tig00004487_arrow/tig00004487_arrow.gff | perl -ne 'my @t=split "\t"; if (@t<8) { next; } elsif (/^#/) { next; } else {print;}' > carAur01.gene.x1.gff
gzip carAur01.gene.x1.gff
zcat  ../maker4/carAur01.gene.gff3.gz  carAur01.gene.x1.gff.gz | perl -e 'while(<>) { if (/^#\s*FASTA\s*$/) { last; } elsif (/^#/) { next; } else { my @t=split /\t/; my $src=$t[1]; if ($src=~m/^\s*$/) {next;} if (!exists $fh{$src}) { open $fh{$src}, "| gzip -c > carAur01.gene.$src.gff.gz"; print {$fh{$src}} "##gff-version 3\n"; } print {$fh{$src}} $_; } } foreach $src (keys(%fh)) {close $fh{$src}; }'
#zcat  carAur01.gene.x1.gff.gz | perl -e 'while(<>) { if (/^#\s*FASTA\s*$/) { last; } elsif (/^#/) { next; } else { my @t=split /\t/; my $src=$t[1]; if ($src=~m/^\s*$/) {next;} if (!exists $fh{$src}) { open $fh{$src}, "| gzip -c > carAur01.gene.$src.x1.gff.gz"; print {$fh{$src}} "##gff-version 3\n"; } print {$fh{$src}} $_; } } foreach $src (keys(%fh)) {close $fh{$src}; }'

for a in maker maker.augustus_masked maker.evm maker.non_overlapping_ab_initio
do
    for f in transcripts proteins
    do
        cat goldfish.arrow.renamed.maker.output/goldfish.arrow.renamed_datastore/B3/74/tig00002911_arrow/tig00002911_arrow.$a.$f.fasta goldfish.arrow.renamed.maker.output/goldfish.arrow.renamed_datastore/81/5B/tig00214208_arrow/tig00214208_arrow.$a.$f.fasta ../maker4/goldfish.arrow.renamed.maker.output/goldfish.arrow.renamed_datastore/FB/57/tig00004487_arrow/tig00004487_arrow.$a.$f.fasta > carAur01.$a.$f.x1.fasta
        cat ../maker4/carAur01.all.$a.$f.fasta carAur01.$a.$f.x1.fasta > carAur01.all.$a.$f.fasta 
        cat carAur01.all.$a.$f.fasta  | sed '/^>/ s/\s.*$//' > carAur01.all.$a.$f.short_name.fasta
    done
done

# mask genome using carAur01.gene.repeatmasker.gff.gz and carAur01.gene.repeatrunner.gff.gz
cat carAur01.gene.repeatmasker.gff carAur01.gene.repeatrunner.gff > carAur01.repeats.merged.gff
bedtools maskfasta -fi ../goldfish.arrow.renamed.fasta -bed carAur01.repeats.merged.gff -soft -fo ../carAur01/carAur01.masked.fasta


##########################
# interproscan
##########################
# interproscan -t p -f TSV,GFF3,HTML -iprlookup --goterms --pathways --minsize 30 carAur01.maker.proteins.x1.fasta interProScan5_run1_x1 400
cat ../maker4/interProScan5_run1/carAur01.ips/*.tsv interProScan5_run1_x1/*.tsv > carAur01.interProScan5.tsv
ips_add_name.pl carAur01.interProScan5.tsv carAur01.interProScan5.add_name.tsv

perl  /home/chenz11/data_cmd/goldfish/11549472/sergey_canu70x/arrow/maker4a/assign_name.pl carAur01.gene.maker.gff.gz t.blast/carAur.ENSDAR.f.join.m6 ~/data/ensembl85/ens85.tid_cdna_ttype_gid_gname_gtype carAur01.gene.maker.add_name.gff
gzip -c carAur01.gene.maker.add_name.gff > carAur01.gene.maker.add_name.gff.gz
perl  ~/data_cmd/goldfish/11549472/sergey_canu70x/arrow/maker4a/combine_maker_ips.pl carAur01.gene.maker.add_name.gff carAur01.interProScan5.add_name.tsv carAur01.gene.maker.add_name_ips.gff
gzip -c carAur01.gene.maker.add_name_ips.gff > carAur01.gene.maker.add_name_ips.gff.gz
echo "##gff-version 3" > ../carAur01/carAur01.gene.gff
gunzip -c carAur01.gene.maker.add_name_ips.gff.gz  >> ../carAur01/carAur01.gene.gff

###########################
# GO stat
###########################
# gene transcript protein  go_domain goterm  goname
cat carAur01.interProScan5.add_name.tsv | awk -F$'\t' '$14~/GO/ {print $1"\t"$14}' | perl -ne '
chomp; @t=split "\t"; @go=split /\|/, $t[1]; $g=$t[0]; $g=~s/-mRNA.*$//;
foreach $go (@go) {@a=split ":",$go; print "$g\t$t[0]\t$t[0]\t$a[2]\t$a[0]:$a[1]\t$a[3]\n"}' |uniq > carAur01.interProScan5.GO.txt
cat carAur01.interProScan5.GO.txt | sort -k4,4 -k5,5 | awk -F$'\t' '{print $4"__"$5}' | uniq -c | awk '{print $2"\t"$1}' | sed 's/__/\t/g' > carAur01.interProScan5.GO.count_by_GO.txt
###########################

cat ~/data/ensembl85/all.85.bed | perl -ne '{chomp; @t=split /\t/; @t1=split /\|/, $t[3]; $t1[0]=substr($t1[0],0,6); if ($t1[2] eq ".") {next;} $t[3]="$t1[0]|$t1[2]"; print join("\t",@t), "\n";}' > sp5.pairs/all.85.bed
cat sp5.pairs/all.85.bed | grep ENSDAR > sp5.pairs/ENSDAR.bed
cat sp5.pairs/all.85.bed | grep ENSAMX > sp5.pairs/ENSAMX.bed
cat ~/data/common_carp/ng/V2.0.gtf | grep '^[^#]' | perl -e '
my $tid0=".", $pid0=".", $gid0=".";
my $chr0, $b0, $e0, $ss0;
my $el0, $cl0;
while(<>) {
    chomp;
    @t = split "\t";
    @attr = split ";",$t[8];
    foreach my $a (@attr) {
        $a=~s/^\s+//; $a=~s/\s+$//;
        my ($u,$v)=split " ", $a;
        $u=~s/^[\s"]+//; $u=~s/[\s"]+$//;
        $v=~s/^[\s"]+//; $v=~s/[\s"]+$//;
        if ($u eq "transcript_id") { $tid=$v; }
        elsif ($u eq "gene_id") { $gid=$v; }
        elsif ($u eq "protein_id") { $pid=$v; }
    }
    if ($t[2] eq "exon") {
        if ($tid ne $tid0) {
        if ($gid0 ne ".") { print "$chr0\t$b0\t$e0\tCYPCAR|$pid0\t1000\t$ss0\n"; }
        $gid0 = $gid;
        $tid0 = $tid;
        $pid0 = $tid;
        $b0   = $t[3]-1;
        $e0   = $t[4];
        $chr0 = $t[0];
        $el0  = 0;
        $cl0  = 0;
        $ss0  = $t[6];
        }
        if ($b0>$t[3]-1) { $b0 = $t[3]-1; }
        if ($e0<$t[4]) { $e0 = $t[4]; }
        $el0 += $t[4]+1-$t[3];
    }
}
if ($gid0 ne ".") { print "$chr0\t$b0\t$e0\tCYPCAR|$pid0\t1000\t$ss0\n"; }
' > sp5.pairs/CYPCAR.bed
cd $cwd

cat ../../carAur03/carAur03.gene.annot.bgp | perl -ne '@t=split "\t"; $t[3]="carAur|$t[3]"; print join("\t",@t[0..5]), "\n";' > sp5.pairs/carAur.bed

zcat ~/data/grass_carp/ng/C_idella_female_genemodels.v1.gmap.gff3.gz | perl -ne '
@t=split "\t"; if ($t[2] ne "gene") {next;} 
@attr = split ";",$t[8];
foreach my $a (@attr) {
    $a=~s/^\s+//; $a=~s/\s+$//;
    my ($u,$v)=split "=", $a;
    $u=~s/^[\s"]+//; $u=~s/[\s"]+$//;
    $v=~s/^[\s"]+//; $v=~s/[\s"]+$//;
    if ($u eq "Name") { $id=$v; }
}
$b = $t[3]-1;
$e = $t[4];
print "$t[0]\t$b\t$e\tCTEIDE|$id\t1000\t$t[6]\n";
' > sp5.pairs/CTEIDE.bed

#--------------------------------------------------------

zcat *.f3.join.m6.gz | awk '$7<$8 && $9<$10 {print $1"\t"$2"\t+\t"$12"\t"$3}' | sed 's/\(ENS....[0-9]\+\)\.[0-9]\+/\1/' | sort -k1,1 -k2,2 | uniq | sort -k1,1 -k4,4gr > all5.f3.join.edge
zcat *.carAur.f4.join.m6.gz | awk '$7<$8 && $9<$10 {print $1"\t"$2"\t+\t"$12"\t"$3}' | sed 's/\(ENS....[0-9]\+\)\.[0-9]\+/\1/' | sort -k1,1 -k2,2 | uniq | sort -k1,1 -k4,4gr > all5.f4.join.edge

cat  ~/data/ensembl85/all.85.bed ~/data/grass_carp/ng/C_idella_female_genemodels.v1.0.LG.bed ~/data/common_carp/ng/V2.0.bed $canu70x_dir/carAur03/carAur03.gene.1.f.bed > $canu70x_dir/carAur03/ens85_GC_CC_GF.bed

