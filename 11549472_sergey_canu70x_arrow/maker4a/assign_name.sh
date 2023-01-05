datadir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/maker4a
cwd=`pwd`

cd $datadir
zcat carAur01.gene.cdna2genome.gff.gz | perl -ne '
    chomp;
	my @tab=split "\t";
	if ($tab[2] eq "expressed_sequence_match") {
		my $name;
		foreach my $attr (split /\s*;\s*/, $tab[8]) {
			my ($u,$v) = split "=", $attr;
			if ($u eq "Name") {
				$name=$v;
			}
		}
		if ($name =~ m/^ENSDART*/) {
			print "$tab[0]\t$tab[3]\t$tab[4]\t$name\t$tab[5]\t$tab[6]\n"
		}
	}
' | gzip -c > carAur01.gene.cdna2genome.danRer10.bed.gz

zcat carAur01.gene.cdna2genome.danRer10.bed.gz | awk -F$'\t' '{ i=index($4,"."); if (i==0) {print $4"\t"$0} else {print substr($4,1,i-1)"\t"$0}}' | sort -k1,1 > carAur01.gene.cdna2genome.danRer10.tmp.a
cat ~/data/ensembl85/ens85.tid_tname_ttype_gid_gname_gtype_zfin_zfinname.txt | cut -f 1,4,5,6 | \
		sort -k1,1 | uniq >  carAur01.gene.cdna2genome.danRer10.tmp.a1
join -t $'\t' -j 1 carAur01.gene.cdna2genome.danRer10.tmp.a carAur01.gene.cdna2genome.danRer10.tmp.a1 \
	| awk -F$'\t' -v OFS=$'\t' '{print $2,$3,$4,$8"__"$9"__"$10,$6,$7}' \ 
	| sort -k1,1 -k2,2n | uniq > carAur01.gene.cdna2genome.danRer10.gene.bed
rm carAur01.gene.cdna2genome.danRer10.tmp.*

zcat carAur01.gene.maker.gff.gz | perl -ne '
    chomp;
	my @tab=split "\t";
	if ($tab[2] eq "gene") {
		my $id;
		foreach my $attr (split /\s*;\s*/, $tab[8]) {
			my ($u,$v) = split "=", $attr;
			if ($u eq "ID") {
				$id=$v;
			}
		}
		print "$tab[0]\t$tab[3]\t$tab[4]\t$id\t$tab[5]\t$tab[6]\n"
	}
'  > carAur01.gene.maker.bed
bedtools intersect -s -f 0.5 -wao -a carAur01.gene.maker.gff.gz \
			 -b carAur01.gene.cdna2genome.danRer10.gene.bed \
		| perl -e 'my @prev; my %name;
		while(<STDIN>) {
			chomp; my @tab=split "\t"; 
			if ($tab[2]!="gene") { print join("\t",@tab[0..8]), "\n"; next; }
			if (@prev<1) { @prev=@tab; next; } 
			if ($tab[0]!=$prev[0] || $tab[3]!=$prev[3] || $tab[4]!=$prev[4] ) {
				if (keys(%name)>0) {
					$prev[8]=~s/([\t;])Name=/\1MakerName=/;
					$prev[8].=";Name=". join(",",keys(%name));
				}
				print join("\t",@prev[0..8]), "\n"; 
				@prev = @tab;
				undef %name;
			}
			if ($tab[12]=~m/[A-Za-z]/) { $name{$tab[12]}++; }
		}
		if (keys(%name)>0) {
			$prev[8]=~s/([\t;])Name=/\1MakerName=/;
			$prev[8].=";Name=". join(",",keys(%name));
		}
		print join("\t",@prev[0..8]), "\n"; ' \
		| sort -k1,6 -k13,13 \


cd $cwd
