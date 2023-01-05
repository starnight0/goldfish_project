cwd=`pwd`
datadir=`echo $cwd | sed 's/data_cmd/data/'`

module load ucsc bedtools

cd $datadir

mkdir bbi;
# hgGcPercent -wigOut -doGaps -file=stdout -win=200 -verbose=0 carAur02 carAur02.2bit | wigEncode stdin gc.wig gc.wib
mv gc.wig gc.wib bbi/
hgGcPercent -wigOut -doGaps -file=stdout -win=200 -verbose=0 carAur02 carAur02.2bit > bbi/gc.wiggle
wigToBigWig bbi/gc.wiggle carAur02.chromSizes bbi/gc.bw


igff=../maker4a/carAur01.gene.maker.add_name_ips.gff
ln -s $igff $asm1.gene.gff
asm1=carAur01
genome1=../goldfish.arrow.renamed.masked
genome1fa=../goldfish.arrow.renamed.masked.fasta
asm2=carAur02
genome2=$asm2
genome2fa=$asm2.fa
chain=carAur01_to_carAur02.liftOver.chain
ogff=$asm.gene.maker.add_name_ips.gff 

echo "##gff-version 3" > $asm1.gene.1.gff 
cat $asm1.gene.gff | perl -ne 'if (m/^#/) { print; next; }
	my @tab = split "\t";
	my $s = ($tab[6] eq "+") ? "F" : "R";
	my $attrs="";
	foreach $attr (split /\s*;\s*/, $tab[8]) {
		my ($u,$v) = split "=", $attr;
		if ($u eq "ID") {
		} elsif ($u eq "Parent") {
		} elsif ($u eq "Name") {
			if ($v =~ m/,\.\.\.$/) { $attr=~s/,\.\.\./\.\.\./; }
			if ($v =~ m/^(ENS.*)__(.*)__(.*)$/) {
				if ($attrs ne "") { $attrs .= ";"; }
				$attrs .= "name2=$2";
			}
		} elsif ($u !~ m/^[a-z]/) {
			$attr="u$attr";
		}
		if ($attrs ne "") { $attrs .= ";"; }
		$attrs .= $attr;
	}
	$tab[8] = $attrs;
	print join("\t",@tab), "\n";
	' >>  $asm1.gene.1.gff 
gff3ToGenePred -attrsOut=$asm1.gene.1.gp.attr $asm1.gene.1.gff $asm1.gene.1.gp
cat $asm1.gene.1.gp.attr | perl -e 'while(<>) {chomp; @t=split "\t"; if ($t[1] eq "name2" && $t[2]=~m/[A-Za-z]/) {if (!exists $a{$t[0]}) {$a{$t[0]}=["",""];} $a{$t[0]}[0]=$t[2];} elsif ($t[1] eq "Name" && $t[2]=~m/[A-Za-z]/) { if (!exists $a{$t[0]}) {$a{$t[0]}=["",""];} $a{$t[0]}[1]=$t[2]} } foreach my $u (keys(%a)) { print "$u\t$a{$u}[0]\t$a{$u}[1]\n"}' > $asm1.gene.1.gp.geneNames
liftOver -multiple -genePred $asm1.gene.1.gp $chain $asm2.gene.1.gp $asm2.gene.1.gp.unmap
genePredToBigGenePred -geneNames=$asm1.gene.1.gp.geneNames $asm2.gene.1.gp stdout | sort -k1,1 -k2,2n > $asm2.gene.1.bgp
bedToBigBed -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as -extraIndex=name,name2,geneName,geneName2 $asm2.gene.1.bgp $asm2.chromSizes $asm2.gene.1.bb

# must load cufflinks
# gffread $asm1.gene.gff -T -o $asm1.gene.gtf
# ldHgGene -genePredExt -out=$asm1.gene.gp -nobin $asm1 gene $asm1.gene.gtf
# liftOver -genePred $asm1.gene.gp $chain $asm2.gene.gp $asm2.gene.gp.unmap
#liftOver -gff $asm1.gene.gtf $chain $asm2.gene.maker.gtf $asm2.gene.maker.gtf.unmap
liftOver -multiple -gff $asm1.gene.gff $chain $asm2.gene.maker.gff $asm2.gene.maker.gff.unmap
#cat $asm2.gene.maker.gtf | perl -ne 'if (m/^#/) { print; next; }
#	my @tab = split "\t";
#	my $s = ($tab[6] eq "+") ? "F" : "R";
#	s/(transcript_id\s*\")([^\-]+)/$1$tab[0]$s-$2/;
#	s/(gene_id\s*\")([^\-]+)/$1$tab[0]$s-$2/;
#	print $_;
#	' >  $asm2.gene.maker.LG.gtf 
#gtfToGenePred -geneNameAsName2 -genePredExt $asm2.gene.maker.LG.gtf $asm2.gene.maker.LG.gp

echo "#gff-version 3" > $asm2.gene.maker.LG.gff 
cat $asm2.gene.maker.gff | perl -ne 'if (m/^#/) { print; next; }
	my @tab = split "\t";
	my $s = ($tab[6] eq "+") ? "F" : "R";
	my $attrs="";
	foreach $attr (split /\s*;\s*/, $tab[8]) {
		my ($u,$v) = split "=", $attr;
		if ($u eq "ID") {
			$attr="ID=$tab[0]$s-$v";
		} elsif ($u eq "Parent") {
			$attr="$u=$tab[0]$s-$v";
		} elsif ($u eq "Name") {
			if ($v =~ m/,\.\.\.$/) { $attr=~s/,\.\.\.//; }
		} elsif ($u !~ m/^[a-z]/) {
			$attr="u$attr";
		}
		if ($attrs ne "") { $attrs .= ";"; }
		$attrs .= $attr;
	}
	$tab[8] = $attrs;
	print join("\t",@tab), "\n";
	' >>  $asm2.gene.maker.LG.gff 
gff3ToGenePred -rnaNameAttr=Name -geneNameAttr=Name $asm2.gene.maker.LG.gff $asm2.gene.maker.LG.gp
cat $asm2.gene.maker.LG.gff | awk '$3=="gene"' | cut -f 9 | perl -ne 'if (m/^#/) {next;}
	chomp;
	my @tab = split ";";
	my $id="";
	my $name="";
	my $name2="";
	for (my $i=0; $i<=$#tab; $i++) {
		if ($tab[$i]=~m/^ID=(.*)$/) {$id=$1;}
		elsif ($tab[$i]=~m/^Name=(.*)$/) {$name=$1;}
		elsif ($tab[$i]=~m/^go=(.*)$/) {$name2=$1;}
	}
	$name2=substr($name2,0,255);
	print "$id\t$name\t\n";
	' > $asm2.gene.maker.gid_name
genePredToBigGenePred -geneNames=$asm2.gene.maker.gid_name $asm2.gene.maker.LG.gp stdout | sort -k1,1 -k2,2n > $asm2.gene.maker.LG.bgp 
genePredToBigGenePred $asm2.gene.maker.LG.gp stdout | sort -k1,1 -k2,2n > $asm2.gene.maker.LG.bgp 
# genePredToBed $asm2.gene.maker.LG.gp $asm2.gene.maker.LG.gp.bed12 
bedToBigBed -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as $asm2.gene.maker.LG.bgp $asm2.chromSizes $asm2.gene.maker.LG.bb

cat $asm2.gene.maker.add_name_ips.gff | awk '$3=="mRNA"' | perl -ne 'if (m/^#/) {next;}
	chomp;
	my @tab = split "\t";
	my @attr=split ";", $tab[8];
	my $id="";
	my $name="";
	my $name2="";
	for (my $i=0; $i<@attr; $i++) {
		if ($attr[$i]=~m/^ID=(.*)$/) {$id=$1;}
		elsif ($attr[$i]=~m/^Name=(.*)$/) {$name=$1;}
		elsif ($attr[$i]=~m/^go=(.*)$/) {$name2=$1;}
	}
	$name2=substr($name2,0,255);
	print "$tab[0]\t", $tab[3]-1, "\t$tab[4]\t$id\t1000\t$tab[6]\n";
	' | sort -k1,1 -k2,2n > $asm2.gene.maker.add_name_ips.bed


ldHgGene -genePredExt -out=$asm2.maker.gene.gp $asm2 gene $asm2.gene.maker.add_name_ips.gff

# cdna2genome track
echo "##gff-version 3" >$asm1.gene.cdna2genome.gff
zcat ../maker4a/carAur01.gene.cdna2genome.gff.gz | perl -ne 'if (m/^#/) { print; next; }
	chomp;
	my @tab = split "\t";
	my $s = ($tab[6] eq "+") ? "F" : "R";
	my %default=(ID=>0, Name=>0, Alias=>0, Parent=>0, Target=>0, Gap=>0, Derives_from=>0, Note=>0, Xref=>0, Ontology_term=>0);
	my $attrs="";
	if ($tab[2] eq "expressed_sequence_match") {$tab[2]="cDNA_match"; }
	foreach $attr (split /\s*;\s*/, $tab[8]) {
		my ($u,$v) = split "=", $attr;
		if (exists $default{$u}) {
		} else {
			$attr="u$attr";
		}
		if ($attrs ne "") { $attrs.=";"; }
		$attrs .= $attr;
	}
	$tab[8] = $attrs;
	print join("\t",@tab), "\n";
	' >> $asm1.gene.cdna2genome.gff
perl ~/my_program3/src/utility/czl_gff_match_part_to_match.pl -i $asm1.gene.cdna2genome.gff -o $asm1.gene.cdna2genome.1.gff
gff3ToPsl $genome1.chromSizes ~/data/ensembl85/ens85_and_carp_gcarp.rna.short_name.chromSizes $asm1.gene.cdna2genome.1.gff $asm1.gene.cdna2genome.psl
pslSwap $asm1.gene.cdna2genome.psl $asm1.gene.cdna2genome.psl.tmp
mv $asm1.gene.cdna2genome.psl.tmp $asm1.gene.cdna2genome.psl
pslPosTarget $asm1.gene.cdna2genome.psl $asm1.gene.cdna2genome.psl.tmp 
mv $asm1.gene.cdna2genome.psl.tmp $asm1.gene.cdna2genome.psl
# add gene name
cat $asm1.gene.cdna2genome.psl | awk '{print $10"\t"$0}' | sort -k1,1 > $asm1.gene.cdna2genome.psl.tmp
sort -k1,1 ~/data/ensembl85/ens85.tid_cdna_ttype_gid_gname_gtype.danRer10 > $asm1.gene.cdna2genome.psl.tmp2
join -t $'\t' -j 1 $asm1.gene.cdna2genome.psl.tmp $asm1.gene.cdna2genome.psl.tmp2 | awk '{a=$2; for (i=3;i<=10; i++) {a=a"\t"$i;} a=a"\t"$11"__"$26; for (i=12;i<=22;i++) {a=a"\t"$i;} print a}' > $asm1.gene.cdna2genome.psl.tmp3 
join -t $'\t' -j 1 -v 1 $asm1.gene.cdna2genome.psl.tmp $asm1.gene.cdna2genome.psl.tmp2 | cut -f 2-22 >> $asm1.gene.cdna2genome.psl.tmp3 
cat $asm1.gene.cdna2genome.psl.tmp3 > $asm1.gene.cdna2genome.psl
rm $asm1.gene.cdna2genome.psl.tmp*
#
pslToBigPsl $asm1.gene.cdna2genome.psl stdout | sort -k1,1 -k2,2n > $asm1.gene.cdna2genome.bigpsl
bedToBigBed -type=bed12+13 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigPsl.as $asm1.gene.cdna2genome.bigpsl $genome1.chromSizes $asm1.gene.cdna2genome.bb
#
liftOver -pslT $asm1.gene.cdna2genome.psl $chain $asm2.gene.cdna2genome.psl $asm2.gene.cdna2genome.psl.unmap
pslToBigPsl $asm2.gene.cdna2genome.psl stdout | sort -k1,1 -k2,2n > $asm2.gene.cdna2genome.bigpsl
bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigPsl.as $asm2.gene.cdna2genome.bigpsl $genome2.chromSizes $asm2.gene.cdna2genome.bb


###########################
# tRNA
fn=$asm2.tRNA
liftOver -bedPlus=12 -tab ../miRBase_run1/carAur01.tRNA.biggp $chain $fn.bgp $fn.bgp.unmap
cat $fn.bgp | sort -k1,1 -k2,2n > $fn.bgp.tmp; 
cat $fn.bgp.tmp | perl -ne '@t=split "\t"; $t[3].="__".$t[12]; print join "\t",@t;' > $fn.bgp
rm $fn.bgp.tmp
bedToBigBed -extraIndex=name2 -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as $fn.bgp $genome2fa.fai $fn.bb
cp $fn.bb ~/data/datashare/fishCarAur1/carAur02/bbi/
###########################

###########################
# miRNA
fn=$asm2.miRNA
liftOver -bedPlus=12 -tab ../miRBase_run1/ex_split/all.out.bigpsl $chain $fn.bigpsl $fn.bigpsl.unmap
cat $fn.bigpsl | sort -k1,1 -k2,2n > $fn.bigpsl.tmp;  mv $fn.bigpsl.tmp $fn.bigpsl
# cat $fn.bigpsl.tmp | perl -ne '@t=split "\t"; $t[3].="__".$t[12]; print join "\t",@t;' > $fn.bigpsl
bedToBigBed -extraIndex=name -type=bed12+13 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigPsl.as $fn.bigpsl $genome2fa.fai $fn.bigpsl.bb
cp $fn.bigpsl.bb ~/data/datashare/fishCarAur1/carAur02/bbi/

###########################

###########################
# miRNA
fn=$asm2.rfam
liftOver -bedPlus=12 -tab ../miRBase_run1/carAur01.rfam.bgp $chain $fn.bgp $fn.bgp.unmap
cat $fn.bgp | sort -k1,1 -k2,2n > $fn.bgp.tmp;  mv $fn.bgp.tmp $fn.bgp
bedToBigBed -extraIndex=name,name2 -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as $fn.bgp $genome2fa.fai $fn.bgp.bb
cp $fn.bgp.bb ~/data/datashare/fishCarAur1/carAur02/bbi/


cd $cwd
