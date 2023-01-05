datadir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/miRBase_run1
cwd=`pwd`

genome1fa=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/goldfish.arrow.renamed.masked.fasta

cd $datadir

cat ex_split/out.*[0-9] > ex_split/all.out
fn=ex_split/all.out
~/my_program3/src/annot_genome/czl_exonerate_to_gtf_gff3.pl -i ex_split/all.out.transcript_exon.gff3 -o ex_split/all.out
cat ex_split/all.out.gff3 | awk '$3=="transcript" || $3=="exon"' > ex_split/all.out.transcript_exon.gff3

~/my_program3/src/annot_genome/czl_exonerate_after_blasn_filter.pl -i ex_split/all.out -o all.ex.out.f -qs ~/data/miRBase/hairpin.U2T.fa.fai -score-frac 0.9 -qcov 0.9 -qiden1 0.8 -qiden2 0.8;
~/my_program3/src/annot_genome/czl_exonerate_to_gtf_gff3.pl -i all.ex.out.f -o all.ex.out.f -qs ~/data/miRBase/hairpin.U2T.fa.fai;
~/my_program3/src/annot_genome/czl_exonerate_after_blasn_gff3_remove_dup.pl -i all.ex.out.f.gff3 -o all.ex.out.f.rmdup.gff;
cuffcompare -o all.ex.f.cuffcompare all.ex.out.f.rmdup.gff 

f=all.ex.f.cuffcompare.combined.gtf
cat $f | perl -e '
while(<>) {
    if (/^#/ || /^\s*$/) {next; } chomp; @t=split "\t"; 
    my $name;
    my $gid;
    foreach $aa (split /\s*;\s*/,$t[8]) {
        my ($u,$v) = split /\s+/, $aa;
        $v=~s/"(.*)"/$1/;
        if ($u eq "gene_id") { $gid=$v; }
        elsif ($u eq "geneID" && !defined $gid) { $gid=$v; }
        if ($u eq "transcript_id") { $tid=$v; }
        if ($u eq "oId") { $v=~s/^[0-9]+_//; $v=~s/^[0-9]+_//; $name=$v; }
    }
    if ($gid=~m/"(.*)"/) { $gid=$1; }
    if ($tid=~m/"(.*)"/) { $tid=$1; }
    print "$gid\t$tid\t$name\t", $_, "\n";
}' > $f.XLOC_TCONS_name

cat $f | perl -e '
while(<>) {
    if (/^#/ || /^\s*$/) {next; } chomp; @t=split "\t"; 
    my $name;
    my $gid;
    foreach $aa (split /\s*;\s*/,$t[8]) {
        my ($u,$v) = split /\s+/, $aa;
        if ($u eq "gene_id") { $gid=$v; }
        elsif ($u eq "geneID" && !defined $gid) { $gid=$v; }
        if ($u eq "transcript_id") { $tid=$v; }
    }
    if ($gid=~m/"(.*)"/) { $gid=$1; }
    if ($tid=~m/"(.*)"/) { $tid=$1; }
    print "$gid\t$tid\t", $_, "\n";
}' > $f.XLOC_TCONS


cat all.ex.f.cuffcompare.combined.gtf | perl -e 'my ($chr0, $b0,$e0,$gid0,$s0);
while(<>) {
    if (/^#/ || /^\s*$/) {next; } chomp; @t=split "\t"; 
    my $name;
    my $gid;
    foreach $aa (split /\s*;\s*/,$t[8]) {
        my ($u,$v) = split /\s+/, $aa;
        if ($u eq "gene_id") { $gid=$v; }
        elsif ($u eq "geneID" && !defined $gid) { $gid=$v; }
    }
    if ($gid=~m/"(.*)"/) { $gid=$1; }
    if (defined $gid0 && $gid0 eq $gid) { if ($b0>$t[3]) {$b0=$t[3];} if ($e0<$t[4])     {$e0=$t[4];} }
    else {
        if (defined $gid0) { print "$chr0\t", $b0-1, "\t$e0\t$gid0\t0\t$s0", "\n"; }
        $chr0=$t[0]; $gid0=$gid; $b0=$t[3]; $e0=$t[4]; $s0=$t[6];
    }
}' > all.ex.f.cuffcompare.combined.gtf.XLOC.bed

bedtools subtract -a "$f.XLOC.bed" -b ../carAur01/carAur01.gene.annot.exon.bed -A | sort -k4,4 > $f.XLOC.not_exon.bed;
join -t$'\t' -1 4 -2 1 "$f.XLOC.not_exon.bed" "$f.XLOC_TCONS" | cut -f 8- > $f.not_exon.gtf;
cat $f.not_exon.gtf | perl -e ' my ($prev_gid, $prev_tid, $prev_name, $gid, $tid, $name); my @ts;
while(<>) {
    if (/^#/ || /^\s*$/) {next; } chomp; my @t=split "\t"; 
    my %info;    my ($chr,$b,$e) = ($t[0], $t[3], $t[4]);
    undef $name;
    foreach $aa (split /\s*;\s*/,$t[8]) {
        my ($u,$v) = split /\s+/, $aa;
        $v=~s/"(.*)"/$1/;
        if ($u eq "gene_id") { $gid=$v; } elsif ($u eq "geneID") { $gid=$v; }
        elsif ($u eq "transcript_id") { $tid=$v; }
        elsif ($u eq "oId") { $v=~s/^[0-9]+_//; $v=~s/^[0-9]+_//; $name=$v; }
        else { $info{$u} = $v; }
    }
    $gid=~s/^XLOC_/CAMIRG/;     $tid=~s/^TCONS_/CAMIRT/;
    if (!defined $prev_gid || $gid ne $prev_gid) {
        if (defined $prev_gid) {
            my @g;
            for (my $i=0; $i<8; $i++) { $g[$i] = $ts[0][$i]; } $g[2] = "gene";
            foreach my $t (@ts) { if ($g[3]>$t->[3]) { $g[3]=$t->[3]; }    if ($g[4]<$t->[4]) { $g[4]=$t->[4]; } }
            print join("\t", @g), "\tID=$prev_gid;Name=$prev_name\n";
            foreach my $t1 (@ts) {
                print join("\t", @{$t1}[0..7]), "\tID=$t1->[8]{ID};Parent=$t1->[8]{Parent};Name=$t1->[8]{Name};exon_num=", $#{$t1->[8]{exons}}+1 , "\n";
                foreach my $e1 (@{$t1->[8]{exons}}) {
                    print join("\t", @$e1[0..7]), "\tID=$e1->[8]{ID};Parent=$e1->[8]{Parent};exon_num=$e1->[8]{exon_num}", "\n";
                }
            }
        }
        $prev_gid = $gid;   $prev_name = $name; @ts=();
    }
    if (!defined $prev_tid || $prev_tid ne $tid) {
        my @t1 = @t[0..7];
        $t1[2]="miRNA";   $t1[8] = {Name=>$name, ID=>$tid, Parent=>$gid, biotype=>"miRNA"};
        push @ts, \@t1;
        $prev_tid = $tid;
    }
    $info{ID}="$tid.E$info{exon_number}";   $info{Parent}=$tid;
    $info{gene_id}=$gid;   $info{exon_num}=$info{exon_number};
    delete $info{exon_number};
    push @{$ts[$#ts][8]{exons}}, [@t[0..7], \%info];
    if ($ts[$#ts][3] > $b) { $ts[$#ts][3] = $b; }     if ($ts[$#ts][4] < $e) { $ts[$#ts][4] = $e; }
}
        if (defined $prev_gid) {
            my @g;
            for (my $i=0; $i<8; $i++) { $g[$i] = $ts[0][$i]; } $g[2] = "gene";
            foreach my $t (@ts) { if ($g[3]>$t->[3]) { $g[3]=$t->[3]; } if ($g[4]<$t->[4]) { $g[4]=$t->[4]; } }
            print join("\t", @g), "\tID=$prev_gid;Name=$prev_name\n";
            foreach my $t1 (@ts) {
                print join("\t", @{$t1}[0..7]), "\tID=$t1->[8]{ID};Parent=$t1->[8]{Parent};Name=$t1->[8]{Name};exon_num=", $#{$t1->[8]{exons}}+1 , "\n";
                foreach my $e1 (@{$t1->[8]{exons}}) {
                    print join("\t", @$e1[0..7]), "\tID=$e1->[8]{ID};Parent=$e1->[8]{Parent};exon_num=$e1->[8]{exon_num}", "\n";
                }
            }
        }
' > $f.not_exon.gff
#gffread -o $f.not_exon.gff -K -M $f.not_exon.gtf

bedtools subtract -a "$f.XLOC.bed" -b ../carAur01/carAur01.gene.annot.expand_1k.bed -A -f 0.5 | sort -k4,4 > $f.XLOC.not_gene.bed;
join -t$'\t' -1 4 -2 1 "$f.XLOC.not_gene.bed" "$f.XLOC_TCONS" | cut -f 8- > $f.not_gene.gtf;
gffread -o $f.not_gene.gff -K -M $f.not_gene.gtf
cat $f.not_gene.gtf | perl -e 'my ($chr0, $b0,$e0,$gid0,$s0);
while(<>) {
    if (/^#/ || /^\s*$/) {next; } chomp; @t=split "\t"; 
    my $name;
    my $gid;
    foreach $aa (split /\s*;\s*/,$t[8]) {
        my ($u,$v) = split /\s+/, $aa;
        if ($v=~m/"(.*)"/) { $v=$1; }
        if ($u eq "gene_id") { $gid=$v; }
        elsif ($u eq "geneID" && !defined $gid) { $gid=$v; }
        elsif ($u eq "oId" && !defined $name) { $name=$v; }
    }
    if ($gid=~m/"(.*)"/) { $gid=$1; }
    if (defined $gid0 && $gid0 eq $gid) { }
    else {
        if ($name=~m/^[0-9]+_[0-9]+_...-(.+)$/i)  {$name=$1;}
        $count{$name}++;
        $gid0 = $gid;
    }
}
foreach my $name (sort keys(%count)) {print "$name\t$count{$name}\n";}' > $f.not_gene.gff.count_by_name.txt

cat $f.not_gene.gtf | perl -e 'my ($chr0, $b0,$e0,$gid0,$s0);
while(<>) {
    if (/^#/ || /^\s*$/) {next; } chomp; @t=split "\t"; 
    my $name;
    my $gid;
    foreach $aa (split /\s*;\s*/,$t[8]) {
        my ($u,$v) = split /\s+/, $aa;
        if ($v=~m/"(.*)"/) { $v=$1; }
        if ($u eq "gene_id") { $gid=$v; }
        elsif ($u eq "geneID" && !defined $gid) { $gid=$v; }
        elsif ($u eq "oId" && !defined $name) { $name=$v; }
    }
    if ($gid=~m/"(.*)"/) { $gid=$1; }
    if (defined $gid0 && $gid0 eq $gid) { }
    else {
        if ($name=~m/^[0-9]+_[0-9]+_(.+)$/i)  {$name=$1;}
        if ($name=~m/dre/) { $count{$name}++; }
        $gid0 = $gid;
    }
}
foreach my $name (sort keys(%count)) {print "$name\t$count{$name}\n";}' > $f.not_gene.gff.count_by_dre_name.txt

# get transcript fasta
gff3ToGenePred ex_split/all.out.transcript_exon.gff3 ex_split/all.out.gp
genePredToBigGenePred ex_split/all.out.gp stdout | sort -k1,1 -k2,2n > ex_split/all.out.bgp
cut -f 1-12 ex_split/all.out.bgp > ex_split/all.out.bed12
bedtools getfasta -s -fi ../goldfish.arrow.renamed.masked.fasta -bed ex_split/all.out.bed12 -tab -split -fo ex_split/all.out.fasta.tab
sort -k1,1 $fn.fasta.tab | uniq | awk '{print ">"$1"\n"$2}' > $fn.fasta

# after run RNAfold
cat ex_split/all.out.rnafold | perl -ne 'chomp;
    if (m/^>(.*)$/) { if (defined $name) {print "$name\t$seq\t$s\n"; }$name=$1; $s=0; $seq="";} 
    elsif ($seq eq "") {$seq=$_}
    elsif (m/([+\-]*[\.0-9]+[0-9]+)/) { if ($1<$s) { $s=$1; } }' > ex_split/all.out.rnafold.tab

cat ex_split/all.out.gff3 | perl -ne 'if (m/^#/) {print; next;}
chomp; @t=split "\t"; @attr=split /\s*;\s*/,$t[8]; 
foreach $a (@attr) { ($u,$v)=split /\s*=\s*/,$a,2; 
    if ($u eq "Parent") {$parent=$v;} 
    elsif ($u eq "Query") {$query=$v;} 
    elsif ($u eq "Align") {@align=split " ", $v;} 
}
if ($t[2] eq "similarity") { next; }
if ($t[2] eq "transcript") { $t[2]="mRNA"; }
if ($t[2] eq "exon") { $t[8]="ID=$parent.$t[3];".$t[8]; }
if ($t[2] eq "intron") { $t[8]="ID=$parent.$t[3];".$t[8]; }
print join("\t",@t), "\n";' > ex_split/all.out.1.gff3

#~/program/EVidenceModeler-1.1.1/EvmUtils/misc/exonerate_gff_to_alignment_gff3.pl ex_split/all.out est > ex_split/all.out.evm.align.gff3
#cat ex_split/all.out.gff3 | perl -ne 'if (m/^#/) {print; next;}
#chomp; @t=split "\t"; @attr=split /\s*;\s*/,$t[8]; 
#foreach $a (@attr) { ($u,$v)=split /\s*=\s*/,$a,2; 
#    if ($u eq "ID") {$id=$v;} 
#    elsif ($u eq "Parent") {$parent=$v;} 
#    elsif ($u eq "Query") {$query=$v;} 
#    elsif ($u eq "Align") {@align=split " ", $v;} 
#}
#if ($t[2] eq "transcript") { $id0=$id; next; }
#if ($t[2] eq "similarity") {
#    $t[2] = "match";
#    $t[8]="ID=$id0;Target=$query $align[1] $align[2]";
#    print join("\t",@t), "\n";
#}' > ex_split/all.out.2.gff3


gff3ToPsl $genome1fa.fai ~/data/miRBase/hairpin.U2T.fa.fai ex_split/all.out.match.gff3 stdout | pslSwap stdin stdout | pslPosTarget stdin ex_split/all.out.psl
# add rnafolder free energy to name field
cat $fn.psl | awk -F$'\t' '{print $14":"$16"-"$17"("$9")""\t"$0}' | sort -k1,1 > $fn.psl.tmp
sort -k1,1 $fn.rnafold.tab > $fn.psl.tmp2
join -t $'\t' -j 1 $fn.psl.tmp $fn.psl.tmp2 | awk '{a=$2; for (i=3;i<=10; i++) {a=a"\t"$i;} a=a"\t"$11"__"$24; for (i=12;i<=22;i++) {a=a"\t"$i;} print a}' > $fn.psl.tmp3 
join -t $'\t' -j 1 -v 1 $fn.psl.tmp $fn.psl.tmp2 | cut -f 2-22 >> $fn.psl.tmp3 
cat $fn.psl.tmp3 > $fn.psl
rm $fn.psl.tmp*
pslToBigPsl $fn.psl stdout | sort -k1,1 -k2,2n > $fn.bigpsl
bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigPsl.as $fn.bigpsl $genome1fa.fai $fn.bb 

cd $cwd
