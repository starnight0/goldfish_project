#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname fileparse);
use czl_io::base_io;
use czl_io::ucsc_chain;
use czl_io::gff;

sub usage()
{
print<<EOF;
Usage:
  $0 <-i INFILE> <-o OUT> <-C 40|50|80> [-g]
Output:
EOF
}

my $wd="$ENV{HOME}/data/goldfish/11549472/sergey_canu70x/arrow";

my $unmasked_gid_file="$wd/carAur01/carAur01.gene.unmasked.gids";
my $masked_gid_file="$wd/carAur01/carAur01.gene.masked.gids";

my $blat_wd="$wd/carAur03/big/blat";
my $gene_cluster_file="$ENV{HOME}/data/goldfish/11549472/sergey_canu70x/arrow/ohnolog/fish4.cluster/czl_ohno_syn.out3.4/rescue_m1.6.cluster.txt";
my $gene_file="$ENV{HOME}/data/goldfish/11549472/sergey_canu70x/arrow/ohnolog/fish4.cluster/czl_ohno_syn.out3.4/anchor";
my $gene_to_cluster_file="$ENV{HOME}/data/goldfish/11549472/sergey_canu70x/arrow/ohnolog/fish4.cluster/czl_ohno_syn.out3.4/rescue_m1.6.anchor_to_cluster.txt";
my $gtp_file="$ENV{HOME}/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/big/validate_lost/all4.gtpnnbbcc";
my $masked_bgp="$wd/carAur03/big/carAur03.noM.gene.masked.bgp";
my $unmasked_bgp="$wd/carAur03/big/carAur03.noM.gene.unmasked.bgp";

#my @sps = qw(ZF GC CC GF);
#my %ncopy = (ZF=>1, GC=>1, CC=>2, GF=>2);
my @sps = qw(ZF GC CC GF);
my %ncopy = (ZF=>1, GC=>1, CC=>2, GF=>2);

my @cluster;
my %gene;
my %rna;
my %sorted_genes;
my %GF_sorted_masked_genes;
my %rG_align; # cDNA to Genome alignment
my $fin;
my $fout;
my $file;
my $in_file;
my $out_file;
my $out_prefix;
my %novo_loci;
my $C='qC80'; # like qC50 or UM300_or_qC50
my $with_gene=0;

if ($#ARGV<0) { usage(); exit 0; }

for (my $k=0; $k<=$#ARGV; $k++) {
	if ($ARGV[$k] =~ m/^(-i)$/) {
		$in_file = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-C)$/) {
		$C = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-o)$/) {
		$out_file = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-g)$/) {
		$with_gene=1;
	} elsif ($ARGV[$k] =~ m/^(-h|--help|-help)$/) {
		usage();
		exit 0;
	} else {
		die "No option '$ARGV[$k]'\n";
	}
}

my $GF_masked_gene_annot = gff::load_gene_annot($masked_bgp, "bgp");
foreach my $g (values( %{$GF_masked_gene_annot->{gene}} )) {
	my $chr = $g->{seqid};
	push @{$GF_sorted_masked_genes{$chr}}, $g;
}
foreach my $chr (keys(%GF_sorted_masked_genes)) {
	@{$GF_sorted_masked_genes{$chr}} = sort {$a->{begin} <=> $b->{begin}} @{$GF_sorted_masked_genes{$chr}};
}


load_gene($gene_file);
#load_cluster($gene_cluster_file);
#load_gtp($gtp_file);
load_rna_to_genome_PSL();
validate_lost();
check_loci_not_overlap();

output();

sub load_gene
{
	$file=shift;
	my $fin=base_io::czl_open($file, "r") or die;
	while(<$fin>) {
		if (m/^#/) { next; }
		if (m/^\s*$/) { next; }
		s/\s*$//;

		my @tab = split /\t/;
		my ($sp, $chr, $begin, $end, $strand, $anchor_idx, $pid, $tid, $gid, $gname, $dup, $pos_idx, $cid, $csz, $cid2, $csz2) = @tab[0..15];
		$gid =~ s/\.[0-9]+$//;
		$tid =~ s/\.[0-9]+$//;
		$pid =~ s/\.[0-9]+$//;
		$rna{$sp}{$tid}{gene_id} = $gid;
		$gene{$sp}{$gid}  = {
			sp     => $sp,
			id     => $gid,
			tid    => $tid,
			pid    => $pid,
			chr    => $chr,
			begin  => $begin ,
			end    => $end   ,
			strand => $strand,
			cid    => $cid,
			csz    => $csz,
			name   => $gname,
		};
		if ($cid>=0) {
			push @{$cluster[$cid]{gene_ids}{$sp}}, $gid;
			$cluster[$cid]{size} = $csz;
		}
	}

# sort by position
	foreach my $sp (@sps) {
		foreach my $g (values( %{$gene{$sp}} )) {
			my $chr = $g->{chr};
			push @{$sorted_genes{$sp}{$chr}}, $g;
		}
		foreach my $chr (keys(%{$sorted_genes{$sp}})) {
			@{$sorted_genes{$sp}{$chr}} = sort {$a->{begin} <=> $b->{begin}} @{$sorted_genes{$sp}{$chr}};
		}
	}
}

sub load_cluster
{
	my $file=shift;
	my $fin=base_io::czl_open($file, "r") or die;
	while(<$fin>) {
		if (m/^#/) { next; }
		if (m/^\s*$/) { next; }
		s/\s*$//;
		my ($sp, $chr, $gid, $gname, $dup, $cid, $csz) = split /\t/, $_;
		$gene{$sp}{$gid}{id}   = $gid;
		$gene{$sp}{$gid}{sp}   = $sp;
		$gene{$sp}{$gid}{chr}  = $chr;
		$gene{$sp}{$gid}{name} = $gname;
		$gene{$sp}{$gid}{cid}  = $cid;
		push @{$cluster[$cid]{gene_ids}{$sp}}, $gid;
		$cluster[$cid]{size} = $csz;
	}
	close $fin;
}

sub load_gtp
{
	my $file=shift;
	my $fin=base_io::czl_open($file, "r");
	while(<$fin>) {
		if (m/^#/) { next; }
		if (m/^\s*$/) { next; }
		s/\s*$//;
		my @tab = split /\t/;
		my $gid = $tab[0]; # gene
		my $tid = $tab[1]; # transcript/RNA
		my $pid = $tab[2]; # protein
		$rna{$tid}{gene_id}=$gid;
	}
	close $fin;
}

sub load_rna_to_genome_PSL
{
	foreach my $sp1 (@sps) {
		foreach my $sp2 (@sps) {
			if ($sp1 eq $sp2) { next; }

	my $file="$blat_wd/${sp2}_genome/${sp1}_cdna/blat.$C.psl";
	my $fin=base_io::czl_open($file, "r") or die;
	while(<$fin>) {
		if (m/^#/) { next; }
		if (m/^\s*$/) { next; }
		s/\s*$//;
		my @tab = split /\t/;
		$tab[ucsc_chain::PSL_QNAME] =~ s/\.[0-9]+$//;
		my $id1 = $tab[ucsc_chain::PSL_QNAME];
		if (!exists $rna{$sp1}{$id1}) { next; }
		my $gid = $rna{$sp1}{$id1}{gene_id};
		push @{ $rG_align{$sp1}{$gid}{$sp2} }, \@tab;
	}
	close $fin;

		}
	}
}

sub validate_lost
{
	my $m=0;
	for (my $cid=0; $cid<@cluster; $cid++) {
		if (!defined $cluster[$cid]) { next; }
		my $c = $cluster[$cid];
		my %sp_ngene;
		my %sp_chr;
		my $n = 0;
		foreach my $sp (@sps) {
			$sp_ngene{$sp} = $#{$c->{gene_ids}{$sp}}+1;
			foreach my $gid (@{$c->{gene_ids}{$sp}}) {
				my $g = $gene{$sp}{$gid};
				$sp_chr{$sp}{$g->{chr}}++;
			}
			if ($sp_ngene{$sp}>0) { $n++; }
		}
		if ($n<2) { next; }
#		if (! exists$sp_ngene{ZF} && ! exists $sp_ngene{GC}) { next; }
#		if ($sp_ngene{ZF}==0 && $sp_ngene{GC}==0) { next; }
		my ($sp0, $sp1, $sp2); # sp0: ZF or GC, #sp1: species with gene missing, #sp2: species with gene number larger than sp1
#		if (defined $sp_ngene{ZF} && $sp_ngene{ZF}==1) { $sp0 = 'ZF'; }
#		elsif (defined $sp_ngene{GC} && $sp_ngene{GC}==1) { $sp0 = 'GC'; }
#		else { next; }

#		my $gid0 = $c->{gene_ids}{$sp0}[0];
#		if ($sp_ngene{CC} < $sp_ngene{GF}) { $sp1='CC'; $sp2='GF'; }
#		else { $sp1='GF'; $sp2='CC'; }

		my %new_sp_ngene = %sp_ngene;
		$c->{new_size} = $c->{size};
		foreach my $sp1 (qw(ZF GC CC GF)) {
			if ($sp_ngene{$sp1}==$ncopy{$sp1}) {next;}
			my @as1;
			foreach my $sp2 (@sps) {
				if ($sp2 eq $sp1) { next; }
				foreach my $gid2 (@{$c->{gene_ids}{$sp2}}) {
					if (! exists $rG_align{$sp2}{$gid2}) {next;}
					if (! exists $rG_align{$sp2}{$gid2}{$sp1}) {next;}
					if (@{$rG_align{$sp2}{$gid2}{$sp1}}==0) {next;}
					my $as = $rG_align{$sp2}{$gid2}{$sp1};
					if (@$as>1) {
						@$as = sort { $a->[ucsc_chain::PSL_MATCH]+$a->[ucsc_chain::PSL_REPMATCH] <=> $b->[ucsc_chain::PSL_MATCH]+$b->[ucsc_chain::PSL_REPMATCH] } @$as;
					}
					for (my $i=0; $i<=$#$as; $i++) {
						my $a = $as->[$i];
						my $chr1 = $a->[ucsc_chain::PSL_TNAME];
						if (exists $sp_chr{$sp1}{$chr1}) { next; } # not on the same chromosome as those in the current cluster
						my $g = find_gene_by_pos($sorted_genes{$sp1}, $a->[ucsc_chain::PSL_TNAME], $a->[ucsc_chain::PSL_TBEGIN], $a->[ucsc_chain::PSL_TEND]);
# locus is not overlap genes in other clusters
# or 
						my $cid1;
						if (defined $g) {
							if (!$with_gene) { next; }
# Let cluster C2 be the merged cluster of current cluster and the overlapped cluster, gene copies (ignore genes of $sp1) of C2 not over 'copies' (ZF:GC:CC:GF=1:1:2:2)
							$cid1 = $g->{cid};
							my $good=1;
							if ($cid1>=0) {
								my $c1 = $cluster[$cid1];
								my %sp_ngene1;
								foreach my $sp (@sps) {
									if ($sp eq $sp1) {next;}
									$sp_ngene1{$sp} = $#{$c1->{gene_ids}{$sp}}+1;
									if ($sp_ngene1{$sp}+$sp_ngene{$sp}>$ncopy{$sp}) {$good=0; last; }
								}
								if ($good) { $m++; }
							}
							if (!$good) {next;}
						}
						if ($sp1 eq "GF" && defined find_gene_by_pos(\%GF_sorted_masked_genes, $a->[ucsc_chain::PSL_TNAME], $a->[ucsc_chain::PSL_TBEGIN], $a->[ucsc_chain::PSL_TEND]) ) {next;}
						if (defined find_novo_loci_by_pos($sp1, $a->[ucsc_chain::PSL_TNAME], $a->[ucsc_chain::PSL_TBEGIN], $a->[ucsc_chain::PSL_TEND]) ) {next;}
						push @as1, [$a, $g];
					}
				}
			}
			if (@as1==0) { next; }

			if (@as1>1) {
				@as1 = reverse sort {$a->[0][ucsc_chain::PSL_MATCH]<=>$b->[0][ucsc_chain::PSL_MATCH]} @as1;
			}
			my $i0=0;
			for (my $i=0; $i<@as1; $i++) {
				my $a = $as1[$i][0];
				my $g = $as1[$i][1];
				my $chr = $a->[ucsc_chain::PSL_TNAME];
				my $beg = $a->[ucsc_chain::PSL_TBEGIN];
				my $end = $a->[ucsc_chain::PSL_TEND];
				my $good=1;

# make sure not on the same chromosome
				if (exists $sp_chr{$sp1}{$chr}) { next; }

# make sure not overlap with other novo loci in the cluster
				my $is_ovl = 0;
				for (my $j=0; $j<$i0; $j++) {
					if (is_tovl($as1[$j][0], $as1[$i][0])) {$is_ovl=1; last;}
				}
				if ($is_ovl) { next; }

				if (!defined $g) {
					if ($new_sp_ngene{$sp1}+1<=$ncopy{$sp1}) {
						push_sort_locus($sp1, $chr, $beg, $end, $cid);
						push @{$c->{novo_loci}{$sp1}}, $a;
						$new_sp_ngene{$sp1}++;
						$sp_chr{$sp1}{$chr}++;
						$c->{new_size}++;
					} else { $good=0; }
				} else {
					if ($g->{cid}<0 || !defined $cluster[$g->{cid}]) {
						if ($new_sp_ngene{$sp1}+1<=$ncopy{$sp1}) {
							push @{$c->{gene_ids}{$sp1}}, $g->{id};
							$new_sp_ngene{$sp1}++;
							$sp_chr{$sp1}{$chr}++;
							$c->{new_size}++;
						} else { $good=0; }
					} else {
						my $cid1 = $g->{cid};
						my $c1 = $cluster[$cid1];
						foreach my $sp (keys(%{$c1->{gene_ids}})) {
							if ($new_sp_ngene{$sp}+$#{$c1->{gene_ids}{$sp}}+1>$ncopy{$sp1}) { $good=0; last; }
						}
						if ($good) {
							foreach my $sp (keys(%{$c1->{gene_ids}})) {
								if ($#{$c1->{gene_ids}{$sp}}<0) { next; }
								push @{$c->{gene_ids}{$sp}}, @{$c->{gene_ids}{$sp}};
								$new_sp_ngene{$sp1}+=$#{$c1->{gene_ids}{$sp}}+1;
								foreach my $gid (@{$c1->{gene_ids}{$sp}}) {
									$gene{$sp}{$gid}{cid}=$cid;
									my $chr1 = $gene{$sp}{$gid}{chr};
									$sp_chr{$sp}{$chr1}++;
									$c->{new_size}++;
								}
							}
							undef $cluster[$cid];
						}
					}
				}
				if ($good) {
					if ($i0!=$i) { $as1[$i0]=$as1[$i]; }
					$i0++;
				}
			}
			$#as1 = $i0-1;
		}

		if ($c->{new_size}>$c->{size}) {
			print "C", $cid+1, " : ";
			foreach my $sp (sort keys(%new_sp_ngene)) {
				if (!exists $sp_ngene{$sp}) {$sp_ngene{$sp}=0; }
			}
			my $i=0;
			foreach my $sp (sort keys(%sp_ngene)) {
				if ($i>0) { print ","; }
				print "$sp:$sp_ngene{$sp}";
				$i++;
			}
			print " ==> ";
			$i=0;
			foreach my $sp (sort keys(%new_sp_ngene)) {
				if ($i>0) { print ","; }
				print "$sp:$new_sp_ngene{$sp}";
				$i++;
			}
			print "\n";
		}

	}
	print "\n$m\n";
}

sub find_gene_by_pos
{
	my ($genes0, $chr, $begin, $end)=@_;
#	my $genes = $sorted_genes{$sp}{$chr};
	my $genes = $genes0->{$chr};
	my $lo = 0;
	my $hi = $#$genes+1;
	my $gi;
	while ($lo<$hi) {
		my $mi = int( ($lo+$hi)/2 );
		if ($genes->[$mi]{end} <= $begin) { $lo = $mi+1; }
		elsif ($genes->[$mi]{begin} >= $end) { $hi = $mi; }
		else { $gi=$mi; last; }
	}
	if (defined $gi) {
		return $genes->[$gi];
	} else {
		return undef;
	}
}

sub find_novo_loci_by_pos
{
	my ($sp, $chr, $begin, $end)=@_;
	my $loci = $novo_loci{$sp}{$chr};
	my $lo = 0;
	my $hi = $#$loci+1;
	my $gi;
	while ($lo<$hi) {
		my $mi = int( ($lo+$hi)/2 );
		if ($loci->[$mi][1] <= $begin) { $lo = $mi+1; }
		elsif ($loci->[$mi][0] >= $end) { $hi = $mi; }
		else { $gi=$mi; last; }
	}
	if (defined $gi) {
		return $gi;
	} else {
		return undef;
	}
}

sub push_sort_locus
{
	my $sp = shift;
	my $chr = shift;
	my $begin = shift;
	my $end = shift;
	my $cid = shift;

	my $loci = $novo_loci{$sp}{$chr};
	my $lo = 0;
	my $hi = $#$loci+1;
	my $gi;
	while ($lo<$hi) {
		my $mi = int( ($lo+$hi)/2 );
		if ($loci->[$mi][0] <= $begin) { $lo = $mi+1; }
		else { $hi = $mi; }
	}
	my $a = [$begin, $end, $cid];
	splice @{$novo_loci{$sp}{$chr}}, $lo, 0, ($a);

	return $a;
}

sub is_tovl
{
	my $a0 = shift; # PSL alignment 1
	my $a1 = shift; # PSL alignment 2
	my $b = $a0->[ucsc_chain::PSL_TBEGIN] > $a1->[ucsc_chain::PSL_TBEGIN] ? $a0->[ucsc_chain::PSL_TBEGIN] : $a1->[ucsc_chain::PSL_TBEGIN];
	my $e = $a0->[ucsc_chain::PSL_TEND  ] < $a1->[ucsc_chain::PSL_TEND  ] ? $a0->[ucsc_chain::PSL_TEND  ] : $a1->[ucsc_chain::PSL_TEND  ];
	if ($b < $e) { return 1; }
	else { return 0; }
}

sub check_loci_not_overlap
{
	foreach my $sp (sort keys(%novo_loci)) {
		foreach my $chr (sort keys(%{$novo_loci{$sp}})) {
			my $a = $novo_loci{$sp}{$chr};
			foreach (my $i=1; $i<=$#$a; $i++) {
				my $a0 = $a->[$i-1];
				my $a1 = $a->[$i];
				if ($a0->[0]>$a1->[0]) { die "Error: novoloci is not sorted.\n"; }
				if ($a0->[1]>$a1->[0]) { die "Error: novoloci is overlap.\n"; }
			}
		}
	}
	foreach my $sp (@sps) {
		foreach my $chr (keys(%{$sorted_genes{$sp}})) {
			my @regions;
			foreach my $g (@{$sorted_genes{$sp}{$chr}}) {
				push @regions, [$g->{begin}, $g->{end}, $g->{id}];
			}
			if (exists $novo_loci{$sp}{$chr}) {
				my $a = $novo_loci{$sp}{$chr};
				foreach (my $i=0; $i<=$#$a; $i++) {
					push @regions, [$a->[$i][0], $a->[$i][1], "$a->[$i][0]-$a->[$i][1]"];
				}
			}
			if (@regions<2) { next; }
			@regions = sort {$a->[0]<=>$b->[0]} @regions;
			for (my $i=1; $i<=$#regions; $i++) {
				my $a0 = $regions[$i-1];
				my $a1 = $regions[$i];
				if ($a0->[1]>$a1->[0] && ($a0->[2]=~m/^[0-9]+-[0-9]+$/ || $a1->[2]=~m/^[0-9]+-[0-9]+$/) ) { die "Error: $sp:$chr:$a1->[2] is overlap with $sp:$chr:$a0->[2].\n"; }
			}
		}
	}
}

sub output
{
	my $fout = base_io::czl_open($out_file."loci.txt", "w");
	my $fout2 = base_io::czl_open($out_file."cluster.txt", "w");
	for (my $cid=0; $cid<@cluster; $cid++) {
		my $c = $cluster[$cid];
		if (!defined $c) { next; }

		my %sp_ngene;
		foreach my $sp (@sps) {
			if (exists $c->{gene_ids}{$sp} && $#{$c->{gene_ids}{$sp}}>=0 ) {
				$sp_ngene{$sp} = $#{$c->{gene_ids}{$sp}}+1;
			}
			if (exists $c->{novo_loci}{$sp} && $#{$c->{novo_loci}{$sp}}>=0 ) {
				$sp_ngene{$sp} += $#{$c->{novo_loci}{$sp}}+1;
			}
		}
		my $sp_size_str="";
		my $csz = 0;
		foreach my $sp (@sps) {
			if (!exists $sp_ngene{$sp}) { next; }
			$csz += $sp_ngene{$sp};
			if ($sp_size_str ne "") { $sp_size_str .= ","; }
			$sp_size_str .= "$sp:$sp_ngene{$sp}";
		}
		print $fout2 "$cid\t$csz\t$sp_size_str\t";
		my @strs;
		foreach my $sp (@sps) {
			if (exists $c->{gene_ids}{$sp}) {
				foreach my $gid (@{$c->{gene_ids}{$sp}}) {
					my $g = $gene{$sp}{$gid};
					print $fout join("\t", ($sp,$g->{chr}, $g->{begin}, $g->{end}, $g->{strand}, $gid, $g->{name}, $cid)), "\n";
					push @strs, "$sp:$g->{chr}:$gid:$g->{name}";
				}
			}
			if (exists $c->{novo_loci}{$sp}) {
				foreach my $a (@{$c->{novo_loci}{$sp}}) {
					my $chr = $a->[ucsc_chain::PSL_TNAME];
					my $beg = $a->[ucsc_chain::PSL_TBEGIN]+1;
					my $end = $a->[ucsc_chain::PSL_TEND];
					my $strand = $a->[ucsc_chain::PSL_STRAND];
					my $qchr = $a->[ucsc_chain::PSL_QNAME];
					my $qbeg = $a->[ucsc_chain::PSL_QBEGIN]+1;
					my $qend = $a->[ucsc_chain::PSL_QEND];
					print $fout join("\t", ($sp,$chr,$beg,$end,$strand, "$chr:$beg:$end", "$qchr:$qbeg:$qend", $cid)), "\n";
					push @strs, "$sp:$chr:$beg-$end:$qchr|$qbeg|$qend";
				}
			}
		}
		print $fout2 join(';', @strs), "\n";
	}
	close $fout;
	close $fout2;
}
