bgenome=broken
genome=$canu70x_dir/goldfish.arrow.renamed.masked

bgenome=$1
map4f=$2
dir=$3
asm_dir=$4
asm=CHRR_genome
if [ "x$3" != "x" ]; then asm_dir=$4; fi  # path
if [ "x$4" != "x" ]; then asm=$5; fi  #  base name

if [ "x$2" == "x" ]
then 
	echo "$0 BrokenGenomePrefix MAP OUTDIR [ASMDIR] [ASM]"
	echo "Like: $0 broken f.map4 ./ /path/to/carAur02 carAur02"
	exit 1
fi

bgenomefa=$bgenome.fasta

if ! [ -d $dir ]; then mkdir $dir; fi
if ! [ -d $asm_dir ]; then mkdir $asm_dir; fi
cout=$dir/chromonomer_nosplit_out
if ! [ -d $cout ]
then
	mkdir $cout
	~/my_program3/src/utility/czl_vcf_liftover.pl -i $map4f -o $dir/all.br.map4 -c $bgenome.fasta.chain
	~/my_program3/src/utility/czl_onemap2_to_Chromonomer_in.pl -io $dir/all.br.map4 -o $dir/chromonomer_in.
	~/program/chromonomer-1.06/chromonomer -p $dir/chromonomer_in.lg -s $dir/chromonomer_in.sam -a $bgenome.fasta.agp -o $dir/chromonomer_nosplit_out --fasta $bgenomefa --disable_splitting > $dir/chromonomer_nosplit.stdout 2>$dir/chromonomer_nosplit.stderr
	mkdir $cout/ctg
	mv $cout/tig* $cout/ctg/

# rename fasta and sequence
	sed '/^>[0-9]/ s/^>/>LG/' $cout/CHRR_integrated.fa > $asm_dir/$asm.fa
	samtools faidx $asm_dir/$asm.fa
	samtools dict -o $asm_dir/$asm.dict $asm_dir/$asm.fa 
fi

#
#tail -n +4 $cout/CHRR_genome.agp | agpToLift -revStrand > $cout/CHRR_genome.agp.liftUp
cat $cout/CHRR_genome.agp | perl -ne 'if (m/^#/) {print $_; next;} my @tab=split /\t/,$_; if ($tab[4] == "W") { if ($tab[5]=~m/^(tig.*arrow)_([0-9]+)_([0-9]+)$/) { $tab[5]=$1; $tab[6]+=$2-1; $tab[7]+=$2-1; print "LG", join "\t",@tab; } else {print "LG", $_;} } else {print "LG", $_;}' > $asm_dir/$asm.LG.agp
sort -k1,1 $asm_dir/$asm.LG.agp > $asm_dir/$asm.tmp.a1
sort -k1,1 $asm_dir/$asm.fa.fai > $asm_dir/$asm.tmp.a2
join -t $'\t' -j 1 $asm_dir/$asm.tmp.a1 $asm_dir/$asm.tmp.a2 | cut -f 1-10 | awk '{print $6"\t"$0}' | sort -k1,1 > $asm_dir/$asm.tmp.a
rm $asm_dir/$asm.tmp.a1 $asm_dir/$asm.tmp.a2
sort -k1,1 $genome.fasta.fai > $asm_dir/$asm.tmp.a3
join -t $'\t' -j 1 $asm_dir/$asm.tmp.a $asm_dir/$asm.tmp.a3 | cut -f 2-12 | sort -k1,1 -k4,4n > $asm_dir/$asm.LG.agp.add_size
rm $asm_dir/$asm.tmp.a3 $asm_dir/$asm.tmp.a

# get not-contained unplaced contig and size
cat $asm_dir/$asm.LG.agp | awk '$5!="U" && $0 !~ /^#/' | cut -f 6 | sort | uniq > $asm_dir/$asm.placed_contig
cat $genome.fasta.fai | cut -f 1,2 | sort > $asm_dir/$asm.tmp.a
join -t $'\t' -v 1 $asm_dir/$asm.tmp.a $asm_dir/$asm.placed_contig > $asm_dir/$asm.unplaced_contig.chromSizes
cut -f 1 $asm_dir/$asm.unplaced_contig.chromSizes > $asm_dir/$asm.unplaced_contig
cat ../goldfish.arrow.renamed.not_contained.fasta.fai | cut -f 1,2 | sort > $asm_dir/$asm.tmp.a
join -t $'\t' -v 1 $asm_dir/$asm.tmp.a $asm_dir/$asm.placed_contig > $asm_dir/$asm.not_contained.unplaced_contig.chromSizes
cp $asm_dir/$asm.LG.agp  $asm_dir/$asm.agp 
cat $asm_dir/$asm.not_contained.unplaced_contig.chromSizes | awk -v OFS=$'\t' '{print $1,"1",$2,"1","W",$1,"1",$2,"+"}' >> $asm_dir/$asm.agp 
cp $asm_dir/$asm.LG.agp.add_size  $asm_dir/$asm.agp.add_size 
cat $asm_dir/$asm.not_contained.unplaced_contig.chromSizes | awk -v OFS=$'\t' '{print $1,"1",$2,"1","W",$1,"1",$2,"+",$2,$2}' >> $asm_dir/$asm.agp.add_size 

cat $asm_dir/$asm.LG.agp.add_size | perl -e 'my $id=0; while(<STDIN>) { if (m/^#/) { next; } 
	chomp;
	my @tab=split "\t";
	my ($qchr, $qbegin, $qend, $idx, $type, $tchr, $tbegin, $tend, $strand, $qsize, $tsize) = @tab[0..10];
	if ($type eq "U") { next; }
	$id++;
	$score = ($tend-$tbegin)*100;
	$tbegin--;
	$qbegin--;
	print "chain $score $tchr $tsize + $tbegin $tend $qchr $qsize";
	if ($strand eq "+") {
		print " + $qbegin $qend"
	} else {
		print " - ", $qsize-$qend, " ", $qsize-$qbegin;
	}
	print " $id\n";
	my $l = $tend-$tbegin;
	print "$l\n\n";
	}
	' > $asm_dir/carAur01_to_$asm.LG.liftOver.chain

cat $asm_dir/$asm.agp.add_size | perl -e 'my $id=0; while(<STDIN>) { if (m/^#/) { next; } 
	chomp;
	my @tab=split "\t";
	my ($qchr, $qbegin, $qend, $idx, $type, $tchr, $tbegin, $tend, $strand, $qsize, $tsize) = @tab[0..10];
	if ($type eq "U") { next; }
	$id++;
	$score = ($tend-$tbegin)*100;
	$tbegin--;
	$qbegin--;
	print "chain $score $tchr $tsize + $tbegin $tend $qchr $qsize";
	if ($strand eq "+") {
		print " + $qbegin $qend"
	} else {
		print " - ", $qsize-$qend, " ", $qsize-$qbegin;
	}
	print " $id\n";
	my $l = $tend-$tbegin;
	print "$l\n\n";
	}
	' > $asm_dir/carAur01_to_$asm.liftOver.chain

# tail -n +4 $asm_dir/$asm.LG.agp | agpToLift -revStrand > $asm_dir/$asm.liftUp
rm $asm_dir/$asm.tmp*

