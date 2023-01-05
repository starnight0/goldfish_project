# after finishing t.blast/blastn.sh
# choose longest RNA or protein for each species
# run:
cat ~/data/danRer10/ensGene.85.20160923.ucsc.represtive.bed | cut -f 4 | cut -d"|" -f 3 | awk '$0!="." {print "ENSDAR|"$0}' > p.blast/ENSDAR.represtive.pid
cat ~/data/ensembl85/Astyanax_mexicanus.AstMex102.85.represtive.bed | cut -f 4 | cut -d"|" -f 3 | awk '$0!="."{print "ENSAMX|"$0}' > p.blast/ENSAMX.represtive.pid

sps10="CTEIDE CYPCAR ENSAMX ENSDAR ENSGAC ENSGMO ENSLAC ENSLOC ENSONI ENSORL ENSPFO ENSTNI ENSTRU ENSXMA carAur"
mybin=~/my_program3/src
for dir1 in t.blast q.blast
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

	dir2=$dir1/sp5.pairs
	if ! [ -d $dir2 ]; then mkdir $dir2; fi
	for sp1 in ENSDAR carAur CYPCAR CTEIDE ENSAMX
	do
	for sp2 in ENSDAR carAur CYPCAR CTEIDE ENSAMX
	do
		if [[ "$sp1" < "$sp2" || "$sp1" == "$sp2" ]]
		then
			echo $sp1.$sp2
			perl $mybin/ohnolog/czl_blast_pair_filter.pl -i $dir1/pairs/$sp1.$sp2.m6.gz -op $dir2/$sp1.$sp2.f --piden-min 30 --ppositive-min 50 --e-max 0.001 --cov-min-min 20 --cov-max-min 75 --bit-frac-min 0.1
			zcat $dir2/$sp1.$sp2.f.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_filter.bit_frac.pl -i - -o - -f 0.2 | sort -t $'\t' -k1,1 -k12,12gr | gzip -c > $dir2/$sp1.$sp2.f2.join.m6.gz
#			zcat $dir2/$sp1.$sp2.f2.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.pl -i - -op $dir2/$sp1.$sp2 --copy $sp1,1:$sp2,2
#			cat $sp1.$sp2.rbh.txt | tail -n +2 | cut -f 9 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}'
#			cat $sp1.$sp2.rbh.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' 'BEGIN {id="";n=0;m=0;x=100;} { if (id!=$1) { if (id!="" && n>1) {print id"\t"n"\t"(x-m);} id=$1; n=0; m=100; x=0;} n=n+1; if ($9<m) {m=$9;} if ($9>x) {x=$9;} } END{ if (n>1) {print id"\t"n"\t"(x-m)}}' > $sp1.$sp2.rbh.diff.txt
#			cat $sp1.$sp2.rbh.txt | tail -n +2 | sort -k1,1 | awk -F$'\t' 'BEGIN {id="";n=0;ids="";} { if (id!=$1) { if (id!="" && n>1) {print ids;} id=$1; n=1; ids=$2;} else { n=n+1; ids=ids"\t"$2} } END{if (n>1) {print ids}}' > $sp1.$sp2.rbh.$sp2.paralog.txt
		fi
	done
	done
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
~/my_program3/src/utility/czl_tab_join.pl -i1 sp5.pairs/carAur.carAur.f.join.m6.gz -i2 carAur01.gene.sp.not_contained.bed -1 1 -2 4 -o - | awk '$1!="." && $21!="."' | cut -f 1-20  > a
~/my_program3/src/utility/czl_tab_join.pl -i1 a -i2 carAur01.gene.sp.not_contained.bed -1 2 -2 4 -o - | awk '$2!="." && $21!="."' | cut -f 1-20  | gzip -c > sp5.pairs/carAur.carAur.f.join.not_contained.m6.gz

# detect heterozygous genes for out carAur
zcat sp5.pairs/carAur.carAur.f.join.not_contained.m6.gz | perl -e '
while(<>) {
	chomp;
	my @t=split "\t";
	if ($t[1] eq $t[0]) {next;}
	my $iden=$t[2];
	my $qcov = $t[18];
	my $tcov = $t[19];
	if ($iden>98 && $qcov>90 && $tcov>90) {print $_, "\n";}
}' | gzip -c > sp5.pairs/carAur.carAur.f.join.poss_hetero.m6.gz

~/my_program3/src/utility/czl_tab_join.pl -i1 sp5.pairs/carAur.carAur.f.join.poss_hetero.m6.gz -i2 carAur01.gene.sp.not_contained.bed -o - -1 1 -2 4 | awk '$1!="." && $21!="."' | cut -f 1-20,21  > a
~/my_program3/src/utility/czl_tab_join.pl -i1 a -i2 ../../carAur01/contig.stat -o - -1 21 -2 2 | cut -f 1-21,24,28 | awk '$1!="." && $22!="."' > a1
~/my_program3/src/utility/czl_tab_join.pl -i1 a1 -i2 carAur01.gene.sp.not_contained.bed -o - -1 2 -2 4 | cut -f 1-23,24 | awk '$1!="." && $24!="."' > b
~/my_program3/src/utility/czl_tab_join.pl -i1 b -i2 ../../carAur01/contig.stat -o - -1 24 -2 2 | cut -f 1-24,27,31 | awk '$1!="." && $25!="."' > b1
cat b1 | awk '$23<0.8 && $26<0.8' | gzip -c > sp5.pairs/carAur.carAur.f.join.poss_hetero.m6.26.gz
zcat sp5.pairs/carAur.carAur.f.join.poss_hetero.m6.26.gz | perl -e '
my %good; my %bad; my %good_chr; my %bad_chr; 
open OUT1, ">sp5.pairs/carAur.hetero.good.tid";
open OUT2, ">sp5.pairs/carAur.hetero.bad.tid";
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
#				$bad_chr{$id2}++;
			}
		}
	} else {
		if (exists $good{$id2}) {
			if (exists $good_chr{$chr1}) { next; }
			elsif (!exists $bad{$id1}) {
				$bad{$id1} = $chr1;
#				$bad_chr{$id1}++;
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
open OUT1, ">sp5.pairs/carAur.hetero.good.chr";
open OUT2, ">sp5.pairs/carAur.hetero.bad.chr";
foreach my $id (sort keys(%good_chr)) { print OUT1 "$id\n"; }
foreach my $id (sort keys(%bad_chr)) { print OUT2 "$id\n"; }
'
~/my_program3/src/utility/czl_tab_join.pl -i1 carAur01.gene.sp.not_contained.bed -i2 sp5.pairs/carAur.hetero.bad.chr -1 1 -2 1 -o - | awk '$1!="." && $7=="."' > carAur01.gene.sp.not_contained.1.bed
# ~/my_program3/src/utility/czl_tab_join.pl -i1 carAur01.gene.sp.not_contained.bed -i2 sp5.pairs/carAur.hetero.bad.tid -1 4 -2 1 -o - | awk '$1!="." && $7=="."' > carAur01.gene.sp.not_contained.1.bed
## get total size and average size of removed heterozygous contigs
~/my_program3/src/utility/czl_tab_join.pl -i1 sp5.pairs/carAur.hetero.bad.chr -i2 ../../carAur01/contig.stat -1 1 -2 2 -o - | awk '$1!="." && $2!="."' | cut -f 1,4 | awk 'BEGIN{l=0;n=0;} {l+=$2; n++} END{print l"\t"l/n}'
~/my_program3/src/utility/czl_tab_join.pl -i1 sp5.pairs/carAur.carAur.f.join.not_contained.m6.gz -i2 carAur01.gene.sp.not_contained.1.bed -o - -1 1 -2 4 | awk '$1!="." && $21!="."' | cut -f 1-20 > a
~/my_program3/src/utility/czl_tab_join.pl -i1 a -i2 carAur01.gene.sp.not_contained.1.bed -o - -1 2 -2 4 | awk '$1!="." && $21!="."' | cut -f 1-20 | gzip -c > sp5.pairs/carAur.carAur.f.join.not_contained.1.m6.gz 

cat /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/carAur03.gene.1.f.bed | cut -f 4 | cut -d"|" -f 3 | awk '$0!="."{print "carAur|"$0}' > p.blast/carAur.represtive.pid

# keep only representive proteins for each gene
dir2=p.blast/sp5.pairs
zcat $dir2/ENSDAR.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSDAR.represtive.pid -1 1 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSDAR.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/ENSDAR.ENSDAR.f3.join.m6.gz 
zcat $dir2/carAur.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 1 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/carAur.carAur.f3.join.m6.gz 
zcat $dir2/ENSAMX.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSAMX.represtive.pid -1 1 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSAMX.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/ENSAMX.ENSAMX.f3.join.m6.gz 
zcat $dir2/ENSDAR.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSDAR.represtive.pid -1 1 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/ENSDAR.carAur.f3.join.m6.gz 
zcat $dir2/ENSAMX.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSAMX.represtive.pid -1 1 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSDAR.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/ENSAMX.ENSDAR.f3.join.m6.gz 
zcat $dir2/ENSAMX.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSAMX.represtive.pid -1 1 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/ENSAMX.carAur.f3.join.m6.gz 

zcat $dir2/CYPCAR.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSDAR.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/CYPCAR.ENSDAR.f3.join.m6.gz 
zcat $dir2/CTEIDE.ENSDAR.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSDAR.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/CTEIDE.ENSDAR.f3.join.m6.gz 
zcat $dir2/CYPCAR.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/CYPCAR.carAur.f3.join.m6.gz 
zcat $dir2/CTEIDE.carAur.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/carAur.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/CTEIDE.carAur.f3.join.m6.gz 
zcat $dir2/CYPCAR.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSAMX.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/CYPCAR.ENSAMX.f3.join.m6.gz 
zcat $dir2/CTEIDE.ENSAMX.f2.join.m6.gz | sed 's/\(ENS\S\+\)\.[0-9]\+\t/\1\t/g' | \
		 ~/my_program3/src/utility/czl_tab_join.pl -i1 - -i2 p.blast/ENSAMX.represtive.pid -1 2 -2 1 -o - | \
		 awk '$1!="." && $21!="."' | cut -f 1-20 | \
		 gzip -c > $dir2/CTEIDE.ENSAMX.f3.join.m6.gz 

cp $dir2/CTEIDE.CYPCAR.f2.join.m6.gz $dir2/CTEIDE.CYPCAR.f3.join.m6.gz 
cp $dir2/CTEIDE.CTEIDE.f2.join.m6.gz $dir2/CTEIDE.CTEIDE.f3.join.m6.gz 
cp $dir2/CYPCAR.CYPCAR.f2.join.m6.gz $dir2/CYPCAR.CYPCAR.f3.join.m6.gz 

# homolog
for dir1 in t.blast q.blast
do
	dir2=$dir1/sp5.pairs
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



# ENSDAR and CYPCAR
sp1=CYPCAR
sp2=ENSDAR
zcat $dir2/$sp1.$sp2.f3.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.pl -i - -op $dir2/rbh/$sp1.$sp2 --copy $sp1,2:$sp2,1
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

# ENSDAR and carAur
sp1=ENSDAR
sp2=carAur
zcat $dir2/$sp1.$sp2.f3.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.pl -i - -op $dir2/rbh/$sp1.$sp2 --copy $sp1,1:$sp2,2
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


# CYPCAR and carAur
sp1=CYPCAR
sp2=carAur
zcat $dir2/$sp1.$sp2.f3.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.pl -i - -op $dir2/rbh/$sp1.$sp2 --copy $sp1,1:$sp2,1
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | cut -f 9 | perl -ne 'chmop; print int($_*10+0.5), "\n";' | \
		sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.iden.hist
cat $dir2/rbh/$sp1.$sp2.rbh.txt | tail -n +2 | perl -ne 'chmop; @t=split /\t/; if ($t[0] eq $t[1]) {next; } if ($t[8]==100) { next; } print int($t[16]*1000/$t[9]+0.5), "\n";' | \
		sort -k1,1gr | uniq -c | awk '{print $2"\t"$1}' > $dir2/rbh/$sp1.$sp2.rbh.positive.hist

	
# CYPCAR and carAur, 2:2
sp1=CYPCAR
sp2=carAur
zcat $dir2/$sp1.$sp2.f3.join.m6.gz | sort -k12,12gr | perl $mybin/ohnolog/czl_blast_pair_rbh.pl -i - -op $dir2/rbh/$sp1.$sp2.22 --copy $sp1,2:$sp2,2
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
zcat sp5.pairs/*.f.join.m6.gz | cut -f 1,2,11 | awk '$1 != $2' > $mcl_dir/seq.abc
mcxload -abc $mcl_dir/seq.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(100)' -o $mcl_dir/seq.mci -write-tab $mcl_dir/seq.tab
mcl $mcl_dir/seq.mci -I 1.4 -odir $mcl_dir
mcl $mcl_dir/seq.mci -I 2 -odir $mcl_dir
mcl $mcl_dir/seq.mci -I 4 -odir $mcl_dir
mcl $mcl_dir/seq.mci -I 6 -odir $mcl_dir
mcxdump -icl $mcl_dir/out.seq.mci.I14 -tabr $mcl_dir/seq.tab -o $mcl_dir/dump.seq.mci.I14
mcxdump -icl $mcl_dir/out.seq.mci.I20 -tabr $mcl_dir/seq.tab -o $mcl_dir/dump.seq.mci.I20
mcxdump -icl $mcl_dir/out.seq.mci.I40 -tabr $mcl_dir/seq.tab -o $mcl_dir/dump.seq.mci.I40
mcxdump -icl $mcl_dir/out.seq.mci.I60 -tabr $mcl_dir/seq.tab -o $mcl_dir/dump.seq.mci.I60


