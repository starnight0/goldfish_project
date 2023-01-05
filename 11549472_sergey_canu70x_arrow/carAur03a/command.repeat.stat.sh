cwd=`pwd`
datadir=`echo $cwd | sed 's/data_cmd/data/'`

module load ucsc bedtools

cd $datadir
mkdir repeat_stat
cd repeat_stat

asm=carAur03
asm_dir=~/data/datashare/fishCarAur1/$asm
mkdir -p $asm_dir/bbi

kent=/data/genome/jksrc_v352/kent/src/hg/lib/
asm1=carAur01
asm1_dir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01
asm2=carAur03
asm2_dir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur03
chain=$asm2_dir/carAur01_to_$asm2.liftOver.chain


repeat_bgp=carAur03.repeatmasker.bgp
cat ../$repeat_bgp | perl -e 'while(<>) {
if (m/^#/ || m/^\s*$/) {next;}
chomp; @t=split "\t", $_, -1;
$chr=$t[0]; $b=$t[1]; $e=$t[2]; $l=$e-$b;
if ($t[12]=~m/species:(.+)\|genus:(.+)\s*$/) {
    $species=$1; $genus=$2;
    if ($genus=~m/^(.+)\/(.+)$/) { $superfamily=$1; $family=$2; }
    else { $superfamily=$genus; $family=$genus; }
    $stat{$superfamily}{$family}{$species}[0]++;
    $stat{$superfamily}{$family}{$species}[1]+=$l;
    $stat1{$superfamily}{$family}[0]++;
    $stat1{$superfamily}{$family}[1]+=$l;
    $stat2{$superfamily}[0]++;
    $stat2{$superfamily}[1]+=$l;
}
}
foreach my $sf (sort keys(%stat)) {
foreach my $fam (sort keys(%{$stat{$sf}})) {
foreach my $sp (sort keys(%{$stat{$sf}{$fam}})) {
    my $a =$stat{$sf}{$fam}{$sp};
    my $a1=$stat1{$sf}{$fam};
    my $a2=$stat2{$sf};
    print $sf,"\t",$fam,"\t",$sp;
    print "\t", $a2->[0], "\t", $a2->[1];
    print "\t", $a1->[0], "\t", $a1->[1];
    print "\t", $a->[0], "\t", $a->[1];
    print "\n";
} } }
' > $repeat_bgp.class_counts.txt;
echo "SuperFamily"$'\t'"Family"$'\t'"Class"$'\t'"SuperFamily Count"$'\t'"SuperFamily bp"$'\t'"Family Count"$'\t'"Family bp"$'\t'"Class Count"$'\t'"Class bp" > $repeat_bgp.class_counts.header.txt;
cut -f 1,4,5 $repeat_bgp.class_counts.txt | uniq > $repeat_bgp.superfamily_counts.txt;
cut -f 1,2,6,7 $repeat_bgp.class_counts.txt | uniq > $repeat_bgp.family_counts.txt;

cat $repeat_bgp.family_counts.txt | awk -F$'\t' -v OFS=$'\t' '{gsub(/-.*$/, "",$2); print $0}' | sort -k1,1 -k2,2 | awk -F$'\t' -v OFS=$'\t' -v a2="" -v n=0 -v m=0 '{if ($2==a2) {n+=$3;m+=$4;} else { if (n>0) {print a1,a2,n,m;} a1=$1;a2=$2; n=$3; m=$4; }} END{print a1,a2,n,m}' > $repeat_bgp.family_counts.compact.txt

# output sequence for each superfamilies
#mkdir repeat_superfamily
#cat carAur03.repeatmasker.bed | grep -v 'Simple' | grep -v 'Low_complexity' | grep -v 'Unknown' | sort -k4,4 -k1,1 -k2,2n > a
#cat a | perl -e '
#my ($c,$sf);
#while(<>) {
#    if (m/^#/ || m/^\s*$/) { next; }
#    chomp;
#    my @t=split "\t",$_;
#    my ($c1, $sf1);
#    if ($t[3]=~m/species:(\S+)\|genus:(\S+)\/(\S+)/) { $c1=$2; $sf1=$3; $f1=$1; }
#    elsif ($t[3]=~m/species:(\S+)\|genus:(\S+)/) { $c1=$2; $sf1="_$c1"; $f1=$1;}
#    else { print STDERR "Not species:genus\n"; next; }
#    if (!defined $c || $c ne $c1 || $sf ne $sf1) {
#        if (defined $c) { close OUT; }
#        $c=$c1; $sf=$sf1;
#        $c=~s/\?/_/g;
#        $sf=~s/[()\?]/_/g;
#        print "$c\t$sf\n";
#        if ( ! -d "repeat_superfamily/$c" ) { mkdir "repeat_superfamily/$c"; }
#        open OUT, ">repeat_superfamily/$c/$sf.bed" or die "Fail to create repeat_superfamily/$c/$sf.bed\n";
#    }
#    print OUT $_ , "\n";
#}
#close OUT;
#';

perl  ~/my_program3/src/annot_genome/czl_repeatmasker_tree.pl -i carAur03.repeatmasker.bed -b ../RM1/RepeatModeler/merged.lib -g carAur03.sm.fa -o repeat_superfamily 2>czl_repeatmasker_tree.stderr

perl  ~/my_program3/src/annot_genome/czl_repeatmasker_seperate_by_libid.pl -i carAur03.repeatmasker.bed -b ../RM1/RepeatModeler/merged.lib -g carAur03.sm.fa -o repeat_id 2>czl_repeatmasker_seperate_by_libid.stderr

cd repeat_superfamily
file1s=r02.clustalw.small.sh
file1l=r02.clustalw.large.sh
file2s=r02.mafft.small.sh
file2l=r02.mafft.large.sh
>$file1s
>$file1l
>$file2s
>$file2l
file3=r03.fasttree.sh
>$file3
for f in `ls */*.merged.fa`
do
    f=`echo $f | sed 's/\.fa$//'`
    n=`cat $f.fa | grep '^>' | wc -l`;
    if [ $n -lt 10 ]; then continue; fi
    if [ $n -lt 300 ]
    then
        echo 'f='$f'; clustalw -INFILE=$f.fa -ALIGN -OUTFILE=$f.clustalw.fa -OUTPUT=FASTA -TYPE=DNA > $f.clustalw.stdout 2> $f.clustalw.stderr' >> $file1s
        echo 'f='$f'; mafft --maxiterate 1000 --genafpair $f.fa > $f.mafft.fa 2> $f.mafft.stderr' >> $file2s;
    else
        echo 'f='$f'; clustalw -INFILE=$f.fa -ALIGN -OUTFILE=$f.clustalw.fa -OUTPUT=FASTA -TYPE=DNA > $f.clustalw.stdout 2> $f.clustalw.stderr' >> $file1l
        echo 'f='$f'; mafft --retree 2 $f.fa > $f.mafft.fa 2> $f.mafft.stderr' >> $file2l
    fi
    echo 'f='$f'; FastTree -gtr -nt $f.clustalw.fa > $f.clustalw.fasttree' >> $file3
    echo 'f='$f'; FastTree -gtr -nt $f.mafft.fa > $f.mafft.fasttree' >> $file3
done
cd ..
# cmd=r02.clustalw.small; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m clustalw -g 8 --partition quick --time=4:00:00 -J $cmd --logdir=$cmd.log -f $cmd.sh
# cmd=r02.clustalw.large; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m clustalw -g 16 --partition norm --time=48:00:00 -J $cmd --logdir=$cmd.log -f $cmd.sh


