# after finishing t.blast/blastn.sh
# choose longest RNA or protein for each species
# run:
# cat ~/data/danRer10/ensGene.85.20160923.ucsc.represtive.bed | cut -f 4 | cut -d"|" -f 3 | awk '$0!="." {print "ENSDAR|"$0}' > p.blast/ENSDAR.represtive.pid
# cat ~/data/danRer10/ensGene.85.20160923.ucsc.represtive.bed | cut -f 4 | cut -d"|" -f 2 | awk '$0!="." {print "ENSDAR|"$0}' > t.blast/ENSDAR.represtive.tid
# cat ~/data/ensembl85/Astyanax_mexicanus.AstMex102.85.represtive.bed | cut -f 4 | cut -d"|" -f 3 | awk '$0!="."{print "ENSAMX|"$0}' > p.blast/ENSAMX.represtive.pid
# cat ~/data/ensembl85/Astyanax_mexicanus.AstMex102.85.represtive.bed | cut -f 4 | cut -d"|" -f 2 | awk '$0!="."{print "ENSAMX|"$0}' > t.blast/ENSAMX.represtive.tid
cat ../carAur01/carAur01.to_remove_tid.txt | awk '{print "carAur|"$1}' > p.blast/carAur01.to_remove_tid.txt 
cat ../carAur01/carAur01.to_remove_tid.txt | awk '{print "carAur|"$1}' > t.blast/carAur01.to_remove_tid.txt 

sps10="CTEIDE CYPCAR ENSAMX ENSDAR ENSGAC ENSGMO ENSLAC ENSLOC ENSONI ENSORL ENSPFO ENSTNI ENSTRU ENSXMA carAur"
mybin=~/my_program3/src
for dir1 in t.blast p.blast
do
    if ! [ -f $dir1/out.m6.gz ]; then  gzip $dir1/out.m6; fi
    if ! [ -d $dir1/pairs ]
    then
        mkdir $dir1/pairs;
        zcat $dir1/out.m6.gz | perl -e 'while(<>) { @t=split "\t"; 
            $q=$t[0]; $t=$t[1]; 
            if ($q=~m/^ENS/) {$qsp=substr $q,0,6;} elsif ($q=~m/^CI/) {$qsp="CTEIDE";} elsif ($q=~m/^CAFS/) {$qsp="CYPCAR";} else {$qsp="carAur";}
            if ($t=~m/^ENS/) {$tsp=substr $t,0,6;} elsif ($t=~m/^CI/) {$tsp="CTEIDE";} elsif ($t=~m/^CAFS/) {$tsp="CYPCAR";} else {$tsp="carAur";}
            $t[0] = "$qsp|$t[0]";
            $t[1] = "$tsp|$t[1]";
            my $sp1 = $qsp;
            my $sp2 = $tsp;
            if ( ($sp1 cmp $sp2) > 0 ) { $sp1=$tsp; $sp2=$qsp; }
            if (!exists $fh{$sp1}{$sp2}) { open $fh{$sp1}{$sp2}, "| gzip -c > '"$dir1"'/pairs/$sp1.$sp2.m6.gz"; $fh{$sp2}{$sp1}=$fh{$sp1}{$sp2}; }
            $t[12]=".";
            print {$fh{$sp1}{$sp2}} join("\t",@t);
        }
        foreach my $sp1 (keys(%fh)) {
            foreach my $sp2 (keys(%{$fh{$sp1}})) {
                close $fh{$sp1}{$sp2};
            }
        }'
    fi

    dir2=$dir1/sp5.pairs.2
    if ! [ -d $dir2 ]
    then
        mkdir $dir2
        for sp1 in ENSDAR carAur CYPCAR CTEIDE ENSAMX
        do
        for sp2 in ENSDAR carAur CYPCAR CTEIDE ENSAMX
        do
            if [[ "$sp1" < "$sp2" || "$sp1" == "$sp2" ]]
            then
                echo $sp1.$sp2
                perl $mybin/ohnolog/czl_blast_pair_filter.pl -i $dir1/pairs/$sp1.$sp2.m6.gz -op $dir2/$sp1.$sp2.f --piden-min 30 --ppositive-min 50 --e-max 0.001 --cov-min-min 50 --cov-max-min 75 --bit-frac-min 0.1
                zcat $dir2/$sp1.$sp2.f.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_filter.bit_frac.pl -i - -o - -f 0.2 | sort -t $'\t' -k1,1 -k12,12gr | gzip -c > $dir2/$sp1.$sp2.f2.join.m6.gz
#            zcat $dir2/$sp1.$sp2.f2.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.pl -i - -op $dir2/$sp1.$sp2 --copy $sp1,1:$sp2,2
#            cat $sp1.$sp2.rbh.txt | tail -n +2 | cut -f 9 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}'
#            cat $sp1.$sp2.rbh.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' 'BEGIN {id="";n=0;m=0;x=100;} { if (id!=$1) { if (id!="" && n>1) {print id"\t"n"\t"(x-m);} id=$1; n=0; m=100; x=0;} n=n+1; if ($9<m) {m=$9;} if ($9>x) {x=$9;} } END{ if (n>1) {print id"\t"n"\t"(x-m)}}' > $sp1.$sp2.rbh.diff.txt
#            cat $sp1.$sp2.rbh.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' 'BEGIN {id="";n=0;ids="";} { if (id!=$1) { if (id!="" && n>1) {print ids;} id=$1; n=1; ids=$2;} else { n=n+1; ids=ids"\t"$2} } END{if (n>1) {print ids}}' > $sp1.$sp2.rbh.$sp2.paralog.txt
            fi
        done
        done
    fi
done

~/my_program3/src/utility/czl_tab_join.pl -i1 ENSDAR.carAur.rbh.carAur.paralog.txt -i2 carAur.carAur.f.join.m6.gz -o aa -1 1,2 -2 1,2
cat aa | awk '$1!="." && $3!="."' | cut -f 3- > ENSDAR.carAur.rbh.carAur.paralog.m8 
cat ENSDAR.carAur.rbh.carAur.paralog.m8 | cut -f 3 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}'

cd t.blast
cat ../carAur01.gene.maker.gff | grep '^[^#]' | perl -ne '
@t=split "\t";
@attr=split ";", $t[8];
if ($t[2] eq "mRNA") {
    foreach $a (@attr) {
        ($u,$v)=split /=/, $a;
        if ($u eq "ID") { $id=$v; }
    }
    $t[3]--;
    print "$t[0]\t$t[3]\t$t[4]\tcarAur|$id\t1000\t$t[6]\n";
}
' | sort -k1,1 -k2,2n > carAur01.gene.sp.bed
~/my_program3/src/utility/czl_tab_join.pl -i1 carAur01.gene.sp.bed -i2 ../../carAur01/contig.not_contained.stat -1 1 -2 2 -o - | awk '$1!="." && $7!="."' | cut -f 1-6 | uniq > carAur01.gene.sp.not_contained.bed
~/my_program3/src/utility/czl_tab_join.pl -i1 sp5.pairs.2/carAur.carAur.f.join.m6.gz -i2 carAur01.gene.sp.not_contained.bed -1 1 -2 4 -o - | awk '$1!="." && $21!="."' | cut -f 1-20  > a
~/my_program3/src/utility/czl_tab_join.pl -i1 a -i2 carAur01.gene.sp.not_contained.bed -1 2 -2 4 -o - | awk '$2!="." && $21!="."' | cut -f 1-20  | gzip -c > sp5.pairs.2/carAur.carAur.f.join.not_contained.m6.gz

# detect heterozygous genes for our carAur
zcat sp5.pairs.2/carAur.carAur.f.join.not_contained.m6.gz | perl -e '
while(<>) {
    chomp;
    my @t=split "\t";
    if ($t[1] eq $t[0]) {next;}
    my $iden=$t[2];
    my $qcov = $t[18];
    my $tcov = $t[19];
    if ($iden>98 && $qcov>90 && $tcov>90) {print $_, "\n";}
}' | gzip -c > sp5.pairs.2/carAur.carAur.f.join.poss_hetero.m6.gz

~/my_program3/src/utility/czl_tab_join.pl -i1 sp5.pairs.2/carAur.carAur.f.join.poss_hetero.m6.gz -i2 carAur01.gene.sp.not_contained.bed -o - -1 1 -2 4 | awk '$1!="." && $21!="."' | cut -f 1-20,21  > a
~/my_program3/src/utility/czl_tab_join.pl -i1 a -i2 ../../carAur01/contig.stat -o - -1 21 -2 2 | cut -f 1-21,24,28 | awk '$1!="." && $22!="."' > a1
~/my_program3/src/utility/czl_tab_join.pl -i1 a1 -i2 carAur01.gene.sp.not_contained.bed -o - -1 2 -2 4 | cut -f 1-23,24 | awk '$1!="." && $24!="."' > b
~/my_program3/src/utility/czl_tab_join.pl -i1 b -i2 ../../carAur01/contig.stat -o - -1 24 -2 2 | cut -f 1-24,27,31 | awk '$1!="." && $25!="."' > b1
cat b1 | awk '$23<0.8 && $26<0.8' | gzip -c > sp5.pairs.2/carAur.carAur.f.join.poss_hetero.m6.26.gz
zcat sp5.pairs.2/carAur.carAur.f.join.poss_hetero.m6.26.gz | perl -e '
my %good; my %bad; my %good_chr; my %bad_chr; 
open OUT1, ">sp5.pairs.2/carAur.hetero.good.tid";
open OUT2, ">sp5.pairs.2/carAur.hetero.bad.tid";
while(<>) {
    chomp;
    @t=split "\t";
    $id1=$t[0]; $id2=$t[1]; $l1 =$t[13]; $l2 =$t[14];
    $chr1 = $t[20]; $chr2 = $t[23];
    if (exists $good{$id1}) {
        if (exists $good{$id2}) {
            next;
        } else {
            if (exists $good_chr{$chr2}) { next; }
            elsif (!exists $bad{$id2}) {
                $bad{$id2} = $chr2;
#                $bad_chr{$id2}++;
            }
        }
    } else {
        if (exists $good{$id2}) {
            if (exists $good_chr{$chr1}) { next; }
            elsif (!exists $bad{$id1}) {
                $bad{$id1} = $chr1;
#                $bad_chr{$id1}++;
            }
        } else {
            if (exists $bad{$id1} || exists $bad_chr{$chr1}) {
                if (exists $bad{$id2} || exists $bad_chr{$chr2}) {
                    next;
                } else {
                    $good{$id2} = $chr2;
                    $good_chr{$chr2}++;
                }
            } else {
                if (exists $bad{$id2} || exists $bad_chr{$chr2}) {
                    $good{$id1} = $chr1;
                    $good_chr{$chr1}++;
                } else {
                    if ($l1>=$l2) {
                        $good{$id1}=$chr1; $good_chr{$chr1}++;
                        $bad{$id2}=$chr2;
                    }
                    else {
                        $good{$id2}=$chr2; $good_chr{$chr2}++;
                        $bad{$id1}=$chr1;
                    }
                }
            }
        }
    }
}
foreach my $id (sort keys(%good)) {
    print OUT1 "$id\t$good{$id}\n";
}
foreach my $id (sort keys(%bad)) {
    my $chr = $bad{$id};
    if (!exists $good_chr{$chr}) {
        print OUT2 "$id\t$bad{$id}\n";
        $bad_chr{$chr}++;
    } else {
        print OUT1 "$id\t$bad{$id}\n";
    }
}
open OUT1, ">sp5.pairs.2/carAur.hetero.good.chr";
open OUT2, ">sp5.pairs.2/carAur.hetero.bad.chr";
foreach my $id (sort keys(%good_chr)) { print OUT1 "$id\n"; }
foreach my $id (sort keys(%bad_chr)) { print OUT2 "$id\n"; }
'
~/my_program3/src/utility/czl_tab_join.pl -i1 carAur01.gene.sp.not_contained.bed -i2 sp5.pairs.2/carAur.hetero.bad.chr -1 1 -2 1 -o - | awk '$1!="." && $7=="."' > carAur01.gene.sp.not_contained.1.bed
# ~/my_program3/src/utility/czl_tab_join.pl -i1 carAur01.gene.sp.not_contained.bed -i2 sp5.pairs.2/carAur.hetero.bad.tid -1 4 -2 1 -o - | awk '$1!="." && $7=="."' > carAur01.gene.sp.not_contained.1.bed
## get total size and average size of removed heterozygous contigs
~/my_program3/src/utility/czl_tab_join.pl -i1 sp5.pairs.2/carAur.hetero.bad.chr -i2 ../../carAur01/contig.stat -1 1 -2 2 -o - | awk '$1!="." && $2!="."' | cut -f 1,4 | awk 'BEGIN{l=0;n=0;} {l+=$2; n++} END{print l"\t"l/n}'
~/my_program3/src/utility/czl_tab_join.pl -i1 sp5.pairs.2/carAur.carAur.f.join.not_contained.m6.gz -i2 carAur01.gene.sp.not_contained.1.bed -o - -1 1 -2 4 | awk '$1!="." && $21!="."' | cut -f 1-20 > a
~/my_program3/src/utility/czl_tab_join.pl -i1 a -i2 carAur01.gene.sp.not_contained.1.bed -o - -1 2 -2 4 | awk '$1!="." && $21!="."' | cut -f 1-20 | gzip -c > sp5.pairs.2/carAur.carAur.f.join.not_contained.1.m6.gz 

cat /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/carAur03.gene.1.f.bed | cut -f 4 | cut -d"|" -f 3 | awk '$0!="."{print "carAur|"$0}' > p.blast/carAur.represtive.pid

# --------------------------------
# filter 3
# keep only representive proteins for each gene
# p.blast
# {{{
dir2=p.blast/sp5.pairs.2
zcat $dir2/ENSDAR.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 1 -2 1 -o -      | awk '$1!="." && $21!="."' | cut -f 1-20     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20       |  gzip -c > $dir2/ENSDAR.ENSDAR.f3.join.m6.gz 
#zcat $dir2/carAur.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 1 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/carAur.carAur.f3.join.m6.gz 
zcat $dir2/ENSAMX.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20      | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 2 -2 1 -o -      | awk '$1!="." && $21!="."' | cut -f 1-20      | gzip -c > $dir2/ENSAMX.ENSAMX.f3.join.m6.gz 

# ENSDAR carAur
zcat $dir2/ENSDAR.carAur.f2.join.m6.gz   |   sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 1 -2 1 -o -    | awk '$1!="." && $21!="."' | cut -f 1-20 > a1
zcat $dir2/ENSDAR.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2
cat a1 a2 | gzip -c > $dir2/ENSDAR.carAur.f3.join.m6.gz 
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/ENSDAR.carAur.f3.join.m6.gz 

 zcat $dir2/ENSAMX.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20      |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 2 -2 1 -o -      |  awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
 zcat $dir2/ENSAMX.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 1 -2 1 -o -      | awk '$1!="." && $21!="."' | cut -f 1-20     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 2 -2 1 -o -       | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
 cat a1 a2 | gzip -c > $dir2/ENSAMX.ENSDAR.f3.join.m6.gz ;

#zcat $dir2/ENSAMX.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSAMX.represtive.pid -1 1 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/ENSAMX.carAur.f3.join.m6.gz 

zcat $dir2/ENSAMX.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20  > a1;
zcat $dir2/ENSAMX.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20  > a2;
cat a1 a2 | gzip -c > $dir2/ENSAMX.carAur.f3.join.m6.gz ;

zcat $dir2/CYPCAR.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
zcat $dir2/CYPCAR.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CYPCAR.ENSDAR.f3.join.m6.gz ;

zcat $dir2/CTEIDE.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
zcat $dir2/CTEIDE.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CTEIDE.ENSDAR.f3.join.m6.gz ;
;
#zcat $dir2/CYPCAR.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/CYPCAR.carAur.f3.join.m6.gz 
#zcat $dir2/CTEIDE.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/CTEIDE.carAur.f3.join.m6.gz 

zcat $dir2/CYPCAR.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
zcat $dir2/CYPCAR.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CYPCAR.ENSAMX.f3.join.m6.gz ;

zcat $dir2/CTEIDE.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
zcat $dir2/CTEIDE.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CTEIDE.ENSAMX.f3.join.m6.gz ;

cp $dir2/CTEIDE.CYPCAR.f2.join.m6.gz $dir2/CTEIDE.CYPCAR.f3.join.m6.gz ;
cp $dir2/CTEIDE.CTEIDE.f2.join.m6.gz $dir2/CTEIDE.CTEIDE.f3.join.m6.gz ;
cp $dir2/CYPCAR.CYPCAR.f2.join.m6.gz $dir2/CYPCAR.CYPCAR.f3.join.m6.gz ;
cp $dir2/CYPCAR.carAur.f2.join.m6.gz $dir2/CYPCAR.carAur.f3.join.m6.gz ;
cp $dir2/CTEIDE.carAur.f2.join.m6.gz $dir2/CTEIDE.carAur.f3.join.m6.gz ;
cp $dir2/carAur.carAur.f2.join.m6.gz $dir2/carAur.carAur.f3.join.m6.gz ;
# }}}
# --------------------------------
# t.blast
# keep only representive proteins for each gene
# {{{
dir2=t.blast/sp5.pairs.2
zcat $dir2/ENSDAR.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 1 -2 1 -o -      | awk '$1!="." && $21!="."' | cut -f 1-20     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20       |  gzip -c > $dir2/ENSDAR.ENSDAR.f3.join.m6.gz 
#zcat $dir2/carAur.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.tid -1 1 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.tid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/carAur.carAur.f3.join.m6.gz 
zcat $dir2/ENSAMX.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20      | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 2 -2 1 -o -      | awk '$1!="." && $21!="."' | cut -f 1-20      | gzip -c > $dir2/ENSAMX.ENSAMX.f3.join.m6.gz 

# ENSDAR carAur
zcat $dir2/ENSDAR.carAur.f2.join.m6.gz   |   sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 1 -2 1 -o -    | awk '$1!="." && $21!="."' | cut -f 1-20 > a1
zcat $dir2/ENSDAR.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2
cat a1 a2 | gzip -c > $dir2/ENSDAR.carAur.f3.join.m6.gz 
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.tid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/ENSDAR.carAur.f3.join.m6.gz 

 zcat $dir2/ENSAMX.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20      |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 2 -2 1 -o -      |  awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
 zcat $dir2/ENSAMX.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 1 -2 1 -o -      | awk '$1!="." && $21!="."' | cut -f 1-20     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 2 -2 1 -o -       | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
 cat a1 a2 | gzip -c > $dir2/ENSAMX.ENSDAR.f3.join.m6.gz ;

#zcat $dir2/ENSAMX.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSAMX.represtive.tid -1 1 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.tid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/ENSAMX.carAur.f3.join.m6.gz 

zcat $dir2/ENSAMX.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20  > a1;
zcat $dir2/ENSAMX.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20  > a2;
cat a1 a2 | gzip -c > $dir2/ENSAMX.carAur.f3.join.m6.gz ;

zcat $dir2/CYPCAR.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
zcat $dir2/CYPCAR.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CYPCAR.ENSDAR.f3.join.m6.gz ;

zcat $dir2/CTEIDE.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
zcat $dir2/CTEIDE.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSDAR.represtive.tid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CTEIDE.ENSDAR.f3.join.m6.gz ;
;
#zcat $dir2/CYPCAR.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.tid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/CYPCAR.carAur.f3.join.m6.gz 
#zcat $dir2/CTEIDE.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
#         ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.tid -1 2 -2 1 -o - | \
#         awk '$1!="." && $21!="."' | cut -f 1-20 | \
#         gzip -c > $dir2/CTEIDE.carAur.f3.join.m6.gz 

zcat $dir2/CYPCAR.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
zcat $dir2/CYPCAR.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CYPCAR.ENSAMX.f3.join.m6.gz ;

zcat $dir2/CTEIDE.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
zcat $dir2/CTEIDE.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'     | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../ENSAMX.represtive.tid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CTEIDE.ENSAMX.f3.join.m6.gz ;

cp $dir2/CTEIDE.CYPCAR.f2.join.m6.gz $dir2/CTEIDE.CYPCAR.f3.join.m6.gz ;
cp $dir2/CTEIDE.CTEIDE.f2.join.m6.gz $dir2/CTEIDE.CTEIDE.f3.join.m6.gz ;
cp $dir2/CYPCAR.CYPCAR.f2.join.m6.gz $dir2/CYPCAR.CYPCAR.f3.join.m6.gz ;
cp $dir2/CYPCAR.carAur.f2.join.m6.gz $dir2/CYPCAR.carAur.f3.join.m6.gz ;
cp $dir2/CTEIDE.carAur.f2.join.m6.gz $dir2/CTEIDE.carAur.f3.join.m6.gz ;
cp $dir2/carAur.carAur.f2.join.m6.gz $dir2/carAur.carAur.f3.join.m6.gz ;
# }}}
# --------------------------------


# --------------------------------
# filter 4
# {{{
#zcat $dir2/carAur.carAur.f3.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'    | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 1 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20     |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 2 -2 1 -o -      |  awk '$1!="." && $21!="."' | cut -f 1-20     | gzip -c > $dir2/carAur.carAur.f4.join.m6.gz ;
#
#zcat $dir2/ENSDAR.carAur.f3.join.m6.gz   |   sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 1 -2 1 -o -    | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
#zcat $dir2/ENSDAR.carAur.f3.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
#cat a1 a2 | gzip -c > $dir2/ENSDAR.carAur.f4.join.m6.gz ;
#
#zcat $dir2/ENSAMX.carAur.f3.join.m6.gz   |   sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 1 -2 1 -o -    | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
#zcat $dir2/ENSAMX.carAur.f3.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
#cat a1 a2 | gzip -c > $dir2/ENSAMX.carAur.f4.join.m6.gz ;
#
#zcat $dir2/CYPCAR.carAur.f3.join.m6.gz   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 1 -2 1 -o -    | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
#zcat $dir2/CYPCAR.carAur.f3.join.m6.gz   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
#cat a1 a2 | gzip -c > $dir2/CYPCAR.carAur.f4.join.m6.gz ;
#
#zcat $dir2/CTEIDE.carAur.f3.join.m6.gz   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 1 -2 1 -o -    | awk '$1!="." && $21!="."' | cut -f 1-20 > a1;
#zcat $dir2/CTEIDE.carAur.f3.join.m6.gz   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur.represtive.pid -1 2 -2 1 -o -     | awk '$1!="." && $21!="."' | cut -f 1-20 > a2;
#cat a1 a2 | gzip -c > $dir2/CTEIDE.carAur.f4.join.m6.gz ;

zcat $dir2/carAur.carAur.f3.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'    | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 1 -2 1 -o -     | awk '$1!="." && $21=="."' | cut -f 1-20     |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 2 -2 1 -o -      |  awk '$1!="." && $21=="."' | cut -f 1-20     | awk '$1!="" && $1!="." && $2!="" && $2!="."' | gzip -c > $dir2/carAur.carAur.f4.join.m6.gz ;

zcat $dir2/ENSDAR.carAur.f3.join.m6.gz   | awk '$1~/carAur/' | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 1 -2 1 -o -    | awk '$1!="." && $21=="."' | cut -f 1-20 > a1;
zcat $dir2/ENSDAR.carAur.f3.join.m6.gz | awk '$2~/carAur/' |sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 2 -2 1 -o -     | awk '$1!="." && $21=="."' | cut -f 1-20 > a2;
cat a1 a2 | $1!="" && $1!="." && $2!="" && $2!="." | gzip -c > $dir2/ENSDAR.carAur.f4.join.m6.gz ;

zcat $dir2/ENSAMX.carAur.f3.join.m6.gz   | awk '$1~/carAur/' |  sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 1 -2 1 -o -    | awk '$1!="." && $21=="."' | cut -f 1-20 > a1;
zcat $dir2/ENSAMX.carAur.f3.join.m6.gz | awk '$1~/carAur/' |sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g'   |    ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 2 -2 1 -o -     | awk '$1!="." && $21=="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/ENSAMX.carAur.f4.join.m6.gz ;

zcat $dir2/CYPCAR.carAur.f3.join.m6.gz   | awk '$1~/carAur/' |   ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 1 -2 1 -o -    | awk '$1!="." && $21=="."' | cut -f 1-20 > a1;
zcat $dir2/CYPCAR.carAur.f3.join.m6.gz   |  awk '$1~/carAur/' |  ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 2 -2 1 -o -     | awk '$1!="." && $21=="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CYPCAR.carAur.f4.join.m6.gz ;

zcat $dir2/CTEIDE.carAur.f3.join.m6.gz   | awk '$1~/carAur/' |   ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 1 -2 1 -o -    | awk '$1!="." && $21=="."' | cut -f 1-20 > a1;
zcat $dir2/CTEIDE.carAur.f3.join.m6.gz   | awk '$1~/carAur/' |   ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../carAur01.to_remove_tid.txt -1 2 -2 1 -o -     | awk '$1!="." && $21=="."' | cut -f 1-20 > a2;
cat a1 a2 | gzip -c > $dir2/CTEIDE.carAur.f4.join.m6.gz ;
cp $dir2/CTEIDE.CTEIDE.f3.join.m6.gz $dir2/CTEIDE.CTEIDE.f4.join.m6.gz ;
cp $dir2/CTEIDE.CYPCAR.f3.join.m6.gz $dir2/CTEIDE.CYPCAR.f4.join.m6.gz ;
cp $dir2/CTEIDE.ENSAMX.f3.join.m6.gz $dir2/CTEIDE.ENSAMX.f4.join.m6.gz ;
cp $dir2/CTEIDE.ENSDAR.f3.join.m6.gz $dir2/CTEIDE.ENSDAR.f4.join.m6.gz ;
cp $dir2/CYPCAR.CYPCAR.f3.join.m6.gz $dir2/CYPCAR.CYPCAR.f4.join.m6.gz ;
cp $dir2/CYPCAR.ENSAMX.f3.join.m6.gz $dir2/CYPCAR.ENSAMX.f4.join.m6.gz ;
cp $dir2/CYPCAR.ENSDAR.f3.join.m6.gz $dir2/CYPCAR.ENSDAR.f4.join.m6.gz ;
cp $dir2/ENSAMX.ENSAMX.f3.join.m6.gz $dir2/ENSAMX.ENSAMX.f4.join.m6.gz ;
cp $dir2/ENSAMX.ENSDAR.f3.join.m6.gz $dir2/ENSAMX.ENSDAR.f4.join.m6.gz ;
cp $dir2/ENSDAR.ENSDAR.f3.join.m6.gz $dir2/ENSDAR.ENSDAR.f4.join.m6.gz ;
# }}}
# --------------------------------

# --------------------------------
# homolog
# --------------------------------
# {{{
for dir1 in t.blast q.blast
do
    dir2=$dir1/sp5.pairs.2
    >$dir2/homolog_counts.txt
    for sp1 in ENSDAR carAur CYPCAR CTEIDE ENSAMX
    do
    for sp2 in ENSDAR carAur CYPCAR CTEIDE ENSAMX
    do
        if [ -f $dir2/$sp1.$sp2.f3.join.m6.gz ]
        then
            echo $sp1 $sp2
            zcat $dir2/$sp1.$sp2.f3.join.m6.gz | cut -f 1 | grep $sp1 | sort | uniq > a
            zcat $dir2/$sp1.$sp2.f3.join.m6.gz | cut -f 2 | grep $sp1 | sort | uniq >> a
            n1=`cat a | sort | uniq | wc -l | cut -f 1`
            if [ "x$sp1" == "x$sp2" ]
            then
                echo $sp1$'\t'$sp2$'\t'$n1$'\t'$n1 >> $dir2/homolog_counts.txt
            else
                zcat $dir2/$sp1.$sp2.f3.join.m6.gz | cut -f 1 | grep $sp2 | sort | uniq > a
                zcat $dir2/$sp1.$sp2.f3.join.m6.gz | cut -f 2 | grep $sp2 | sort | uniq >> a
                n2=`cat a | sort | uniq | wc -l | cut -f 1`
                echo $sp1$'\t'$sp2$'\t'$n1$'\t'$n2 >> $dir2/homolog_counts.txt
            fi
            zcat $dir2/$sp1.$sp2.f3.join.m6.gz | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[2]==100) { next; } print int($t[17]*1000/$t[3]+0.5), "\n";' | \
                    sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/$sp1.$sp2.positive.hist
            cat $dir2/$sp1.$sp2.positive.hist | perl -e '@h=(0)x1001;
                while(<>) {@t=split "\t"; $h[$t[0]]=$t[1]}
                $w=2;
                @h1=(0)*1001;
                my $i_max=0;
                my $h_max=0;
                for ($i=0;$i<=1000;$i++) {
                    $h1[$i]=0;
                    $n=0;
                    for ($j=$i-$w; $j<=$i+$w;$j++) {
                        if ($j>=0 && $j<=1000) { $h1[$i]+=$h[$j]; $n++;}
                    }
                    $h1[$i]/=$n;
                    if ($h1[$i]>$h_max) {$i_max=$i; $h_max=$h1[$i];}
                }
                print $i_max, "\n" '
        fi
    done
    done
done
# }}}
# --------------------------------


# ENSDAR and CYPCAR
sp1=CYPCAR
sp2=ENSDAR
zcat $dir2/$sp1.$sp2.f4.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.v1.pl -i - -op $dir2/rbh/$sp1.$sp2.3I --copy $sp1,2:$sp2,1 -a $dir2/../all5.6.bed -m 3I
zcat $dir2/$sp1.$sp2.f4.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.v1.pl -i - -op $dir2/rbh/$sp1.$sp2.2I --copy $sp1,2:$sp2,1 -a $dir2/../all5.6.bed -m 2I
zcat $dir2/$sp1.$sp2.f4.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.v1.pl -i - -op $dir2/rbh/$sp1.$sp2.1I --copy $sp1,2:$sp2,1 -a $dir2/../all5.6.bed -m 1I
tail -n +2 $dir2/rbh/$sp1.$sp2.1.rbh.txt | awk '$1~/'$sp2'/' | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 1 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$(NF-1); for (i=1;i<=NF-2;i++) {a=a"\t"$i;} print a}' | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 3 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$1"\t"$2"\t"$(NF-1); for (i=3;i<=NF-2;i++) {a=a"\t"$i;} print a}'  > $dir2/rbh/$sp1.$sp2.1.rbh.with_chr.txt 
tail -n +2 $dir2/rbh/$sp1.$sp2.2.rbh.txt | awk '$1~/'$sp2'/' | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 1 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$(NF-1); for (i=1;i<=NF-2;i++) {a=a"\t"$i;} print a}' | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 3 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$1"\t"$2"\t"$(NF-1); for (i=3;i<=NF-2;i++) {a=a"\t"$i;} print a}'  > $dir2/rbh/$sp1.$sp2.2.rbh.with_chr.txt 
tail -n +2 $dir2/rbh/$sp1.$sp2.3.rbh.txt | awk '$2~/'$sp2'/' | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 1 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$(NF-1); for (i=1;i<=NF-2;i++) {a=a"\t"$i;} print a}' | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 3 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$1"\t"$2"\t"$(NF-1); for (i=3;i<=NF-2;i++) {a=a"\t"$i;} print a}'  | sort -k4,4 > $dir2/rbh/$sp1.$sp2.3.rbh.with_chr.txt

cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | cut -f 9 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.iden.hist
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[8]==100) { next; } print int($t[16]*1000/$t[9]+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.positive.hist
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | sort -k2,2 | awk -F$'\t' 'BEGIN {id="";n=0;m=0;x=100;} { 
        if (id!=$2) { 
            if (id!="" && n>1) {print id"\t"n"\t"(x-m);} 
            id=$2; n=0; m=100; x=0;
        } 
        n=n+1; if ($9<m) {m=$9;} if ($9>x) {x=$9;} } 
        END{ if (n>1) {print id"\t"n"\t"(x-m)}}' > $dir2/rbh/$sp1.$sp2.rbh.diff.txt
cat $dir2/rbh/$sp1.$sp2.rbh.diff.txt | cut -f 3 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.diff.hist.txt
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | sort -k2,2 | awk -F$'\t' 'BEGIN {id="";n=0;ids="";} { 
        if (id!=$2) { 
            if (id!="" && n>1) {print ids;} 
            id=$2; n=1; ids=$1;
        } else { n=n+1; ids=ids"\t"$1 } 
    } END{if (n>1) {print ids}}' > $dir2/rbh/$sp1.$sp2.rbh.$sp1.paralog.txt
~/my_program3/src/utility/czl_tab_join.pl -i1 $dir2/rbh/$sp1.$sp2.rbh.$sp1.paralog.txt -i2 $dir2/$sp1.$sp1.f.join.m6.gz -o aa -1 1,2 -2 1,2
cat aa | awk '$1!="." && $3!="."' | cut -f 3- > $dir2/rbh/$sp1.$sp2.rbh.$sp1.paralog.m6 
cat $dir2/rbh/$sp1.$sp2.rbh.$sp1.paralog.m6 | cut -f 3 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.$sp1.paralog.iden.hist.txt
cat $dir2/rbh/$sp1.$sp2.rbh.$sp1.paralog.m6 | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[2]==100) { next; } $p=int($t[17]*1000/$t[3]+0.5); print $p, "\n";' | \
    sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.$sp1.paralog.positive.hist.txt

cat rbh/$sp1.$sp2.1.rbh.with_chr.txt | perl -ne 'chomp; @t=split "\t"; if ($t[0]!~m/^[0-9]/) {next;} if ($t[2]=~m/^LG([0-9]+)$/) {$c1=$1;} else {next;} if ( int(($c1+1)/2)==$t[0]) {print $_,"\n"}' | sort -k2,2 > rbh/CYPCAR.1.syn2
cat rbh/$sp1.$sp2.2.rbh.with_chr.txt | perl -ne 'chomp; @t=split "\t"; if ($t[0]!~m/^[0-9]/) {next;} if ($t[2]=~m/^LG([0-9]+)$/) {$c1=$1;} else {next;} if ( int(($c1+1)/2)==$t[0]) {print $_,"\n"}' | sort -k2,2 > rbh/CYPCAR.2.syn2
cat rbh/$sp1.$sp2.3.rbh.with_chr.txt | perl -ne 'chomp; @t=split "\t"; if ($t[2]!~m/^[0-9]/) {next;} if ($t[0]=~m/^LG([0-9]+)$/) {$c1=$1;} else {next;} if ( int(($c1+1)/2)==$t[2]) {print $_,"\n"}' | sort -k4,4 > rbh/CYPCAR.3.syn2
cat rbh/CYPCAR.3.syn2 | awk -F$'\t' 'BEGIN{d=0; a1; a2; id="";chr2="";n=0;m=0;} {if (id==$4){n++; if (chr2!=$1) {a2=$0; m++;} } else { if (n==2 && m==2) {print a1;print a2} id=$4; chr2=$1; n=1; m=1; a1=$0;}} END{if (n==2 && m==2) {print a1;print a2}}' > 4bh/CYPCAR.3.syn2.2

# ENSDAR and carAur
sp1=ENSDAR
sp2=carAur
zcat $dir2/$sp1.$sp2.f4.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.v1.pl -i - -op $dir2/rbh/$sp1.$sp2 --copy $sp1,1:$sp2,2 -a $dir2/../all5.6.bed -m 3I
~/my_program3/src/utility/czl_tab_join.pl -i1 $dir2/rbh/$sp1.$sp2.rbh.txt -i2 $dir2/../id_to_chr -1 1 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$(NF-1); for (i=1;i<=NF-2;i++) {a=a"\t"$i;} print a}' | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 3 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$1"\t"$2"\t"$(NF-1); for (i=3;i<=NF-2;i++) {a=a"\t"$i;} print a}'  > $dir2/rbh/$sp1.$sp2.rbh.with_chr.txt 
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | cut -f 9 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.iden.hist
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[8]==100) { next; } print int($t[16]*1000/$t[9]+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.positive.hist
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' 'BEGIN {id="";n=0;m=0;x=100;} { 
        if (id!=$1) { 
            if (id!="" && n>1) {print id"\t"n"\t"(x-m);} 
            id=$1; n=0; m=100; x=0;
        } 
        n=n+1; if ($9<m) {m=$9;} if ($9>x) {x=$9;} } 
        END{ if (n>1) {print id"\t"n"\t"(x-m)}}' > $dir2/rbh/$sp1.$sp2.rbh.diff.txt
cat $dir2/rbh/$sp1.$sp2.rbh.diff.txt | cut -f 3 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.diff.hist.txt
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' 'BEGIN {id="";n=0;ids="";} { 
        if (id!=$1) { 
            if (id!="" && n>1) {print ids;} 
            id=$1; n=1; ids=$2;
        } else { n=n+1; ids=ids"\t"$2 } 
    } END{if (n>1) {print ids}}' > $dir2/rbh/$sp1.$sp2.rbh.$sp2.paralog.txt
~/my_program3/src/utility/czl_tab_join.pl -i1 $dir2/rbh/$sp1.$sp2.rbh.$sp2.paralog.txt -i2 $dir2/$sp2.$sp2.f3.join.m6.gz -o - -1 1,2 -2 1,2 | awk '$1!="." && $3!="."' | cut -f 3- > $dir2/rbh/$sp1.$sp2.rbh.$sp2.paralog.m6 
cat $dir2/rbh/$sp1.$sp2.rbh.$sp2.paralog.m6 | cut -f 3 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.$sp2.paralog.iden.hist.txt
cat $dir2/rbh/$sp1.$sp2.rbh.$sp2.paralog.m6 | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[2]==100) { next; } $p=int($t[17]*1000/$t[3]+0.5); print $p, "\n";' | \
    sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.$sp2.paralog.positive.hist.txt

cat rbh/$sp1.$sp2.rbh.with_chr.txt | perl -ne 'chomp; @t=split "\t"; if ($t[0]!~m/^[0-9]/) {next;} if ($t[2]=~m/^LG([0-9]+)/) {$c1=$1;} else {next;} if ($c1==$t[0] || $c1==$t[0]+25) {print $_,"\n"}' | sort -k2,2 > rbh/carAur.syn2
cat rbh/carAur.syn2 | awk -F$'\t' 'BEGIN{d=0; a1; a2; id="";chr2="";n=0;m=0;} {if (id==$2){n++; if (chr2!=$3) {a2=$0; m++;} } else { if (n==2 && m==2) {print a1;print a2} id=$2; chr2=$3; n=1; m=1; a1=$0;}} END{if (n==2 && m==2) {print a1;print a2}}' > rbh/carAur.syn2.2


# CYPCAR and carAur
sp1=CYPCAR
sp2=carAur
zcat $dir2/$sp1.$sp2.f4.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.v1.pl -i - -op $dir2/rbh/$sp1.$sp2.3 --copy $sp1,1:$sp2,1 -a $dir2/../all5.6.bed -m 3B;
zcat $dir2/$sp1.$sp2.f4.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.v1.pl -i - -op $dir2/rbh/$sp1.$sp2.1B --copy $sp1,1:$sp2,1 -a $dir2/../all5.6.bed -m 1B;
tail -n +2 $dir2/rbh/$sp1.$sp2.rbh.txt | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 1 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$(NF-1); for (i=1;i<=NF-2;i++) {a=a"\t"$i;} print a}' | ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 $dir2/../id_to_chr -1 3 -2 2 -o - | awk -F$'\t' '$1!="." && $NF!="." {a=$1"\t"$2"\t"$(NF-1); for (i=3;i<=NF-2;i++) {a=a"\t"$i;} print a}'  > $dir2/rbh/$sp1.$sp2.rbh.with_chr.txt ;
cat $dir2/rbh/$sp1.$sp2.rbh.with_chr.txt | perl -e 'my %count; my %chr1h; my %chr2h; while(<>) {
    chomp;
    if (m/^#/) {next; }
    my @t=split "\t",$_,-1;
    my $chr1;
    my $chr2;
    if ($t[0]=~m/LG([0-9]+)/) {
        $chr1=$1;
        $chr1h{$chr1}++;
    }
    if ($t[2]=~m/LG([0-9]+)/) {
        $chr2=$1;
        $chr2h{$chr2}++;
    }
    if (defined $chr1 && defined $chr2) {
        $count{$chr1}{$chr2}++;
    }
}
my @chr1s = sort {$a<=>$b} keys(%chr1h);
my @chr2s = sort {$a<=>$b} keys(%chr2h);
my @chr1o = @chr1s;
my @chr2o = @chr2s;
foreach my $chr (@chr1o) { $chr = "LG$chr"; }
foreach my $chr (@chr2o) { $chr = "LG$chr"; }
print join("\t", @chr2o), "\n";
foreach my $chr1 (sort {$a<=>$b} keys(%count)) {
    print "LG$chr1";
    foreach my $chr2 (@chr2s) {
        if (!exists $count{$chr1}{$chr2}) { $count{$chr1}{$chr2}=0; }
        print "\t", $count{$chr1}{$chr2};
    }
    print "\n";
}
' > $dir2/rbh/$sp1.$sp2.rbh.chr_count.txt 
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | cut -f 9 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.iden.hist
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[8]==100) { next; } print int($t[16]*1000/$t[9]+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.positive.hist

# match CC and GF chromosome
cat $dir2/rbh/$sp1.$sp2.3I.rbh.txt | perl -e 'while(<>) {
chomp; @t=split "\t", $_, -1; $chr1=$t[2]; $chr2=$t[3];
if ($chr1=~m/^LG([0-9]+)/) { $i1=$1; } else {next;}
if ($chr2=~m/^LG([0-9]+)/) { $i2=$1; } else {next;}
$count[$i1][$i2]{n}++;
$count[$i1][$i2]{total_iden}+=$t[10];
$count[$i1][$i2]{total_bit} +=$t[15];
}
open OUT1, ">'$dir2'/rbh/'$sp1'.'$sp2'.3I.rbh.pair_count_matrix.txt";
open OUT2, ">'$dir2'/rbh/'$sp1'.'$sp2'.3I.rbh.pair_avg_iden_matrix.txt";
open OUT3, ">'$dir2'/rbh/'$sp1'.'$sp2'.3I.rbh.pair_avg_bit_matrix.txt";
open OUT4, ">'$dir2'/rbh/'$sp1'.'$sp2'.3I.rbh.pair_total_bit_matrix.txt";
open OUT5, ">'$dir2'/rbh/'$sp1'.'$sp2'.3I.rbh.pair_match.txt";
print OUT1 "\t", join("\t", 1..50), "\n"; print OUT2 "\t", join("\t", 1..50), "\n";
print OUT3 "\t", join("\t", 1..50), "\n"; print OUT4 "\t", join("\t", 1..50), "\n";
my %count;
for (my $i1=1; $i1<=50; $i1++) {
    print OUT1 $i1; print OUT2 $i1; print OUT3 $i1; print OUT4 $i1;
	$i2_max=0;
	$i2_max_bit=0;
	$i2_max_iden=0;
    for (my $i2=1; $i2<=50; $i2++) {
        my $aa=$count[$i1][$i2]; 
        if (!defined $aa) {
            print OUT1 "\t0"; print OUT2 "\t0"; print OUT3 "\t0"; print OUT4 "\t0";
            print "$i1\t$i2\t0\t0\t0\t0\n";
        } else {
            $n=$aa->{n}; $iden=$aa->{total_iden}/$n; $bit=$aa->{total_bit}/$n;
			if ($aa->{total_bit}>$i2_max_bit) { $i2_max=$i2; $i2_max_bit=$aa->{total_bit}; $i2_max_iden=$iden; }
            print "$i1\t$i2\t$n\t$iden\t$aa->{total_bit}\t$bit\n";
            print OUT1 "\t$n"; print OUT2 "\t$iden"; print OUT4 "\t$aa->{total_bit}"; print OUT3 "\t$bit";
        }
    }
    print OUT1 "\n"; print OUT2 "\n"; print OUT3 "\n"; print OUT4 "\n";
}
for (my $j=1; $j<=25; $j++) {
	my $i11=$j*2-1; my $i12=$j*2; my $i21=$j; my $i22=$j+25;
	my $c1 = $count[$i11][$i21]{total_bit} + $count[$i12][$22]{total_bit};
	my $c2 = $count[$i11][$i22]{total_bit} + $count[$i12][$21]{total_bit};
	print STDERR "$j\t$i11\t$i12\t$i21\t$i22\t$c1\t$c2\n";
	if ($c1  >= $c2) {
		my $aa=$count[$i11][$i21];
		$n=$aa->{n}; $iden=$aa->{total_iden}/$n; $bit=$aa->{total_bit}/$n;
		print OUT5 "$i11\t$i21\t$n\t$iden\t$aa->{total_bit}\t$bit\n";
		$aa=$count[$i12][$i22];
		$n=$aa->{n}; $iden=$aa->{total_iden}/$n; $bit=$aa->{total_bit}/$n;
		print OUT5 "$i12\t$i22\t$n\t$iden\t$aa->{total_bit}\t$bit\n";
	} else {
		my $aa=$count[$i11][$i22];
		$n=$aa->{n}; $iden=$aa->{total_iden}/$n; $bit=$aa->{total_bit}/$n;
		print OUT5 "$i11\t$i22\t$n\t$iden\t$aa->{total_bit}\t$bit\n";
		$aa=$count[$i12][$i21];
		$n=$aa->{n}; $iden=$aa->{total_iden}/$n; $bit=$aa->{total_bit}/$n;
		print OUT5 "$i12\t$i21\t$n\t$iden\t$aa->{total_bit}\t$bit\n";
	}
}
close OUT1; close OUT2; close OUT3; close OUT4; close OUT5;
' > $dir2/rbh/$sp1.$sp2.3I.rbh.pair_count.txt
    
# CYPCAR and carAur, 2:2
sp1=CYPCAR
sp2=carAur
zcat $dir2/$sp1.$sp2.f4.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.pl -i - -op $dir2/rbh/$sp1.$sp2.22 --copy $sp1,2:$sp2,2
cat $dir2/rbh/$sp1.$sp2.22.rbh.txt | tail -n +2 | cut -f 9 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.22.rbh.iden.hist
cat $dir2/rbh/$sp1.$sp2.22.rbh.txt | tail -n +2 | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[8]==100) { next; } print int($t[16]*1000/$t[9]+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.22.rbh.positive.hist

# CYPCAR and carAur
sp1=CTEIDE
sp2=ENSDAR
zcat $dir2/$sp1.$sp2.f3.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.pl -i - -op $dir2/rbh/$sp1.$sp2 --copy $sp1,1:$sp2,1
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | cut -f 9 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.iden.hist
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[8]==100) { next; } print int($t[16]*1000/$t[9]+0.5), "\n";' | \
        sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.positive.hist

# combine histogram
perl -e 'for (my $i=1000; $i>=700; $i--) {print "$i\n";}' > a
~/my_program3/src/utility/czl_tab_join.pl -i1 a -i2 $dir2/CYPCAR.carAur.rbh.positive.hist -1 1 -2 1 -o - | awk '$1!="." && $NF!="."' | cut -f 1,3- | sort -k1,1nr > a1
~/my_program3/src/utility/czl_tab_join.pl -i1 a1 -i2 $dir2/CYPCAR.ENSDAR.rbh.CYPCAR.paralog.positive.hist.txt -1 1 -2 1 -o - | awk '$1!="." && $NF!="."' | cut -f 1-2,4 | sort -k1,1nr > a2
~/my_program3/src/utility/czl_tab_join.pl -i1 a2 -i2 $dir2/ENSDAR.carAur.rbh.carAur.paralog.positive.hist.txt -1 1 -2 1 -o - | awk '$1!="." && $NF!="."' | cut -f 1-3,5 | sort -k1,1nr > a3
~/my_program3/src/utility/czl_tab_join.pl -i1 a3 -i2 $dir2/CTEIDE.CYPCAR.rbh.positive.hist -1 1 -2 1 -o - | awk '$1!="." && $NF!="."' | cut -f 1-4,6 > a4
~/my_program3/src/utility/czl_tab_join.pl -i1 a4 -i2 $dir2/CTEIDE.carAur.rbh.positive.hist -1 1 -2 1 -o - | awk '$1!="." && $NF!="."' | cut -f 1-5,7 > a5
~/my_program3/src/utility/czl_tab_join.pl -i1 a5 -i2 $dir2/CYPCAR.ENSDAR.rbh.positive.hist -1 1 -2 1 -o - | awk '$1!="." && $NF!="."' | cut -f 1-6,8 > a6
~/my_program3/src/utility/czl_tab_join.pl -i1 a6 -i2 $dir2/ENSDAR.carAur.rbh.positive.hist -1 1 -2 1 -o - | awk '$1!="." && $NF!="."' | cut -f 1-7,9 > a7
~/my_program3/src/utility/czl_tab_join.pl -i1 a7 -i2 $dir2/CTEIDE.ENSDAR.rbh.positive.hist -1 1 -2 1 -o - | awk '$1!="." && $NF!="."' | cut -f 1-8,10 > a8
echo CC_GF$'\t'CC_CC$'\t'GF_GF$'\t'CC_GC$'\t'GF_GC$'\t'CC_ZF$'\t'GF_ZF$'\t'GC_ZF > p.blast/rbh/positive.hist.txt
cat a8 | sort -k1,1gr >> p.blast/rbh/positive.hist.txt

perl -e 'for (my $i=1000; $i>=1; $i--) {print "$i\n";}' > a
~/my_program3/src/utility/czl_tab_join.pl -i1 a -i2 $dir2/CYPCAR.carAur.rbh.iden.hist -1 1 -2 1 -e 0 -o - |  cut -f 1,3- | sort -k1,1gr > a1
~/my_program3/src/utility/czl_tab_join.pl -i1 a1 -i2 $dir2/CYPCAR.ENSDAR.rbh.CYPCAR.paralog.iden.hist.txt -e 0 -1 1 -2 1 -e 0 -o - | cut -f 1-2,4 | sort -k1,1gr > a2
~/my_program3/src/utility/czl_tab_join.pl -i1 a2 -i2 $dir2/ENSDAR.carAur.rbh.carAur.paralog.iden.hist.txt -e 0 -1 1 -2 1 -o - |  cut -f 1-3,5 | sort -k1,1gr > a3
~/my_program3/src/utility/czl_tab_join.pl -i1 a3 -i2 $dir2/CTEIDE.CYPCAR.rbh.iden.hist -e 0 -1 1 -2 1 -o - | cut -f 1-4,6 | sort -k1,1gr > a4
~/my_program3/src/utility/czl_tab_join.pl -i1 a4 -i2 $dir2/CTEIDE.carAur.rbh.iden.hist -e 0 -1 1 -2 1 -o - | cut -f 1-5,7 | sort -k1,1gr > a5
~/my_program3/src/utility/czl_tab_join.pl -i1 a5 -i2 $dir2/CYPCAR.ENSDAR.rbh.iden.hist -e 0 -1 1 -2 1 -o - | cut -f 1-6,8 | sort -k1,1gr > a6
~/my_program3/src/utility/czl_tab_join.pl -i1 a6 -i2 $dir2/ENSDAR.carAur.rbh.iden.hist -e 0 -1 1 -2 1 -o - | cut -f 1-7,9 | sort -k1,1gr > a7
~/my_program3/src/utility/czl_tab_join.pl -i1 a7 -i2 $dir2/CTEIDE.ENSDAR.rbh.iden.hist -e 0 -1 1 -2 1 -o - | cut -f 1-8,10 | sort -k1,1gr > a8
echo Identity$'\t'CC_GF$'\t'CC_CC$'\t'GF_GF$'\t'CC_GC$'\t'GF_GC$'\t'CC_ZF$'\t'GF_ZF$'\t'GC_ZF > p.blast/iden.hist.txt
cat a8 | sort -k1,1gr >> p.blast/iden.hist.txt

######################################################
# MCL cluster 5 species: zebrafish, grass carp, common carp, goldfish, cave fish
######################################################
mcl_dir=sp5.mcl.noself
mkdir $mcl_dir
# Do NOT include self-alignment
zcat sp5.pairs.2/*.f.join.m6.gz | cut -f 1,2,11 | awk '$1 != $2' > $mcl_dir/seq.abc
mcxload -abc $mcl_dir/seq.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(100)' -o $mcl_dir/seq.mci -write-tab $mcl_dir/seq.tab
mcl $mcl_dir/seq.mci -I 1.4 -odir $mcl_dir
mcl $mcl_dir/seq.mci -I 2 -odir $mcl_dir
mcl $mcl_dir/seq.mci -I 4 -odir $mcl_dir
mcl $mcl_dir/seq.mci -I 6 -odir $mcl_dir
mcxdump -icl $mcl_dir/out.seq.mci.I14 -tabr $mcl_dir/seq.tab -o $mcl_dir/dump.seq.mci.I14
mcxdump -icl $mcl_dir/out.seq.mci.I20 -tabr $mcl_dir/seq.tab -o $mcl_dir/dump.seq.mci.I20
mcxdump -icl $mcl_dir/out.seq.mci.I40 -tabr $mcl_dir/seq.tab -o $mcl_dir/dump.seq.mci.I40
mcxdump -icl $mcl_dir/out.seq.mci.I60 -tabr $mcl_dir/seq.tab -o $mcl_dir/dump.seq.mci.I60


