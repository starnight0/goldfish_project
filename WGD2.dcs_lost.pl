#!/usr/bin/perl
#use threads;
#use threads::shared;
use strict;
use Set::IntervalTree;
use warnings;
use Data::Dumper;
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname fileparse);
use czl_io::base_io;
use czl_io::align2;
use czl_io::ucsc_chain;
use czl_io::gff;

sub usage()
{
print<<EOF;
Usage:
  $0 <-i1 GF_ZF.net.chain> <-i2 GF_GF.net.chain> <-o OUT>
Output:
EOF
}

my $can_use_threads = eval 'use threads; 1';

my $fin;
my $fout;
my $file;
my $in_file;
my $in_file1a;
my $in_file1b;
my $in_file2;
my $out_file;
my $out_prefix;
my ($i,$j,$k);
my $exon_file1;
my $exon_file2;
my $nt=8;

my ($tname, $tsize, $tbegin, $tend, $tstrand);
my ($qname, $qsize, $qbegin, $qend, $qstrand);

use constant {
	LOST_BOTH => 0,
	LOST_ONE  => 1,
	LOST_TWO  => 2,
};

if ($#ARGV<0) { usage(); exit 0; }

for (my $k=0; $k<=$#ARGV; $k++) {
	if ($ARGV[$k] =~ m/^(-i1a)$/) {
		$in_file1a = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-i1b)$/) {
		$in_file1b = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-i2)$/) {
		$in_file2 = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-exon1)$/) {
		$exon_file1 = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-exon2)$/) {
		$exon_file2 = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-o)$/) {
		$out_prefix = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-h|--help|-help)$/) {
		usage();
		exit 0;
	} else {
		die "No option '$ARGV[$k]'\n";
	}
}

my @threads;
my @aligns;
my @trees;
my @na;
print STDERR "Load Chains\n";
$aligns[0] = ucsc_chain::load_all_align_from_chain($in_file1a, 0);
$aligns[1] = ucsc_chain::load_all_align_from_chain($in_file1b, 0);
$aligns[2] = ucsc_chain::load_all_align_from_chain($in_file2, 0);

print STDERR "Build Interval Tree\n";
for (my $k=0; $k<3; $k++) {
	$na[$k] = @{$aligns[$k]};
	for (my $i=0; $i<2; $i++) {
		$trees[$k][$i] = align2::build_interval_tree_for_liftover_m2($aligns[$k], $i);
	}
}

#my $v3 = $trees[2][0]{LG26}->fetch(23757279, 23757823);
my @align_matches;

my @exons;
$exons[0] = gff::czl_load_bed_chr_array_hash($exon_file1);
$exons[1] = gff::czl_load_bed_chr_array_hash($exon_file2);
# build interval tree for exons/CNEs
my @exon_trees;
for (my $k=0; $k<2; $k++) {
foreach my $name (keys(%{$exons[$k]})) {
	my $tree = Set::IntervalTree->new;
for (my $ei=0; $ei<@{$exons[$k]{$name}}; $ei++) {
	my $e = $exons[$k]{$name}[$ei];
	$tree->insert([$ei], $e->[1], $e->[2]);
}
	$exon_trees[$k]{$name} = $tree;
}
}
#

print STDERR "Liftover exons\n";
my @map_exons;
#foreach my $name (keys(%{$exons[0]})) {
#	$map_exons[0]{$name} = align2::liftover_m2($exons[0]{$name}, $aligns[0], 1, $trees[0][1]);
#	$map_exons[1]{$name} = align2::liftover_m2($exons[0]{$name}, $aligns[1], 1, $trees[1][1]);
#}
$map_exons[0] = align2::liftover_m3($exons[0], $aligns[0], 1);
$map_exons[1] = align2::liftover_m3($exons[0], $aligns[1], 1);
#foreach my $name (keys(%{$exons[1]})) {
#	$map_exons[2]{$name} = align2::liftover_m2($exons[1]{$name}, $aligns[0], 0, $trees[0][0]);
#	$map_exons[3]{$name} = align2::liftover_m2($exons[1]{$name}, $aligns[1], 0, $trees[1][0]);
#}
$map_exons[2] = align2::liftover_m3($exons[1], $aligns[0], 0);
$map_exons[3] = align2::liftover_m3($exons[1], $aligns[1], 0);
#foreach my $name (keys(%{$exons[1]})) {
#	$map_exons[4]{$name} = align2::liftover_m2($exons[1]{$name}, $aligns[2], 0, $trees[2][0]);
#	$map_exons[5]{$name} = align2::liftover_m2($exons[1]{$name}, $aligns[2], 1, $trees[2][1]);
#}
$map_exons[4] = align2::liftover_m3($exons[1], $aligns[2], 0);
$map_exons[5] = align2::liftover_m3($exons[1], $aligns[2], 1);

print STDERR "Search Tripletes\n";
for (my $k=0; $k<6; $k++) {
foreach my $name (keys(%{$map_exons[$k]})) {
	my $map_exons1 = $map_exons[$k]{$name};
	for (my $ei=0; $ei<@$map_exons1; $ei++) {
		if (@{$map_exons1->[$ei]}==0) { next; }
		if (@{$map_exons1->[$ei]}==1) {
			my $a1 = $map_exons1->[$ei][0];
			$a1->[7] = $a1->[3]-$a1->[2];
			next;
		}
		my @u=(0);
		my $i;
		for ($i=1; $i<@{$map_exons1->[$ei]}; $i++) {
			my $a1 = $map_exons1->[$ei][$i-1];
			my $a2 = $map_exons1->[$ei][$i];
			if ($a2->[0]!=$a1->[0]) { push @u,$i; }
		}
		my @es;
		push @u, $i;
		for (my $ui=1; $ui<@u; $ui++) {
			my ($i0,$i1) = ($u[$ui-1], $u[$ui]);
			my @e = @{$map_exons1->[$ei][$i0]};
			$e[7] = 0;
			for (my $i=$i0; $i<$i1; $i++) {
				my $e1 = $map_exons1->[$ei][$i];
				if ($e1->[2]<$e[2]) { $e[2] = $e1->[2]; }
				if ($e1->[3]>$e[3]) { $e[3] = $e1->[3]; }
				if ($e1->[4]<$e[4]) { $e[4] = $e1->[4]; }
				if ($e1->[5]>$e[5]) { $e[5] = $e1->[5]; }
				$e[7] += $e1->[3]-$e1->[2];
			}
			push @es, \@e; # e: [chain_idx, first_match_in_chain, original_begin, orig_end, mapped_begin, mapped_end, match_length]
		}
		undef $map_exons1->[$ei];
		$map_exons1->[$ei] = \@es;
	}
}
}

foreach my $name (keys(%{$exons[0]})) {
	my $exons1 = $exons[0]{$name};
for (my $i1=0; $i1<@{$exons1}; $i1++) {
	my $exon1 = $exons1->[$i1];
	my @map_exons1 = ($map_exons[0]{$name}[$i1], $map_exons[1]{$name}[$i1]);
	for (my $k1=0; $k1<2; $k1++) {
		for (my $i=0; $i<@{$map_exons1[$k1]}; $i++) {
			my $a1 = $map_exons1[$k1][$i];
			my $aa1 = $aligns[$k1][$a1->[0]];
			my $name2 = align2::get_tname($aa1);
			$a1->[8] = [];
			if (exists $exon_trees[1]{$name2}) {
				my $u = $exon_trees[1]{$name2}->fetch($a1->[4], $a1->[5]);
				if (@$u>0) {
					foreach my $u1 (@$u) {
						push @{$a1->[8]}, $u1->[0];
					}
				}
			}
		}
	}
}
}

foreach my $name (keys(%{$exons[1]})) {
	my $exons1 = $exons[1]{$name};
for (my $i1=0; $i1<@{$exons1}; $i1++) {
	my $exon1 = $exons1->[$i1];
	my @map_exons1 = ($map_exons[2]{$name}[$i1], $map_exons[3]{$name}[$i1]);
	for (my $k1=0; $k1<2; $k1++) {
		for (my $i=0; $i<@{$map_exons1[$k1]}; $i++) {
			my $a1 = $map_exons1[$k1][$i];
			my $aa1 = $aligns[$k1][$a1->[0]];
			my $name2 = align2::get_qname($aa1);
			$a1->[8] = [];
			if (exists $exon_trees[0]{$name2}) {
				my $u = $exon_trees[0]{$name2}->fetch($a1->[4], $a1->[5]);
				if (@$u>0) {
					foreach my $u1 (@$u) {
						push @{$a1->[8]}, $u1->[0];
					}
				}
			}
		}
	}
}
}

foreach my $name (keys(%{$exons[1]})) {
	my $exons1 = $exons[1]{$name};
for (my $i1=0; $i1<@{$exons1}; $i1++) {
	my $exon1 = $exons1->[$i1];
	my @map_exons1 = ($map_exons[4]{$name}[$i1], $map_exons[5]{$name}[$i1]);
	for (my $k1=0; $k1<2; $k1++) {
		for (my $i=0; $i<@{$map_exons1[$k1]}; $i++) {
			my $a1 = $map_exons1[$k1][$i];
			my $aa1 = $aligns[2][$a1->[0]];
			my $name2;
			if ($k1==0) { $name2 = align2::get_qname($aa1); }
			else { $name2 = align2::get_tname($aa1); }
			$a1->[8] = [];
			if (exists $exon_trees[1]{$name2}) {
				my $u = $exon_trees[1]{$name2}->fetch($a1->[4], $a1->[5]);
				if (@$u>0) {
					foreach my $u1 (@$u) {
						push @{$a1->[8]}, $u1->[0];
					}
				}
			}
		}
	}
}
}

# find triple
my @match_loci_count = (0) x 8;
my @exon_stat;
foreach my $name (keys(%{$exons[0]})) {
	my $exons1 = $exons[0]{$name};
for (my $i1=0; $i1<@{$exons1}; $i1++) {
	my $exon1 = $exons1->[$i1];
	my $map_exons1 = $map_exons[0]{$name}[$i1];
	my $map_exons2 = $map_exons[1]{$name}[$i1];
	if (@$map_exons1==0 && @$map_exons2==0) {
		$match_loci_count[0]++;
		$exon_stat[0]{$name}[$i1]{copy_type} = 0;
		next;
	}
	if (@$map_exons1<1 || @$map_exons2<1) {
		$match_loci_count[1]++;
		$exon_stat[0]{$name}[$i1]{copy_type} = 1;
		next;
	}
	if (@$map_exons1>=2 || @$map_exons2>=2) {
		$match_loci_count[3]++;
		$exon_stat[0]{$name}[$i1]{copy_type} = 3;
		next;
	}
	my $j1=0;
	my $j2=0;
	my $a1 = $map_exons1->[$j1];
	my $a2 = $map_exons2->[$j2];
	my $ai1 = $a1->[0];
	my $ai2 = $a2->[0];
	my $aa1 = $aligns[0][$ai1];
	my $aa2 = $aligns[1][$ai2];
	my $name1 = align2::get_tname($aa1);
	my $name2 = align2::get_tname($aa2);
	if ($name1 eq $name2) {
		$match_loci_count[6]++;
		$exon_stat[0]{$name}[$i1]{copy_type} = 6;
		next;
	}
#	my $reg1 = [ $a1->[4], $a1->[5] ];
#	my $reg2 = [ $a2->[4], $a2->[5] ];
# map reg1 to reg3 (the other copy) use alignment3: sp2--sp2 self alifnment
# check overlap of reg3 and reg2
	my $v = [];
	if (exists $trees[2][0]{$name1}) {
		$v = $trees[2][0]{$name1}->fetch($a1->[4], $a1->[5]);
	}
	if (@$v==0) {
		if (exists $trees[2][1]{$name1}) {
			$v = $trees[2][1]{$name1}->fetch($a1->[4], $a1->[5]);
		}
	}
	my ($ai3, $i3, $tbeg3, $tend3, $qbeg3, $qend3);
	if (@$v==0) {
		$match_loci_count[4]++;
		$exon_stat[0]{$name}[$i1]{copy_type} = 4;
		next;
	} elsif (@$v==1) {
		($ai3, $i3, $tbeg3, $tend3, $qbeg3, $qend3) = @{$v->[0]};
	} else {
		my %by_ai;
		foreach my $v1 (@$v) {
			my ($ai4, $i4, $tbeg4, $tend4, $qbeg4, $qend4) = @{$v1};
			if (!exists $by_ai{$ai4}) {
				$by_ai{$ai4} = [$ai4, $i4, $tbeg4, $tend4, $qbeg4, $qend4];
			} else {
				if ($by_ai{$ai4}[2]>$tbeg4) {$by_ai{$ai4}[2]=$tbeg4;}
				if ($by_ai{$ai4}[3]<$tend4) {$by_ai{$ai4}[3]=$tend4;}
				if ($by_ai{$ai4}[4]>$qbeg4) {$by_ai{$ai4}[4]=$qbeg4;}
				if ($by_ai{$ai4}[5]<$qend4) {$by_ai{$ai4}[5]=$qend4;}
			}
		}
		if (keys(%by_ai)>=2) {
			$match_loci_count[7]++;
			$exon_stat[0]{$name}[$i1]{copy_type} = 7;
			next;
		}
		($ai3) = keys(%by_ai);
		($i3,$tbeg3,$tend3,$qbeg3,$qend3) = @{$by_ai{$ai3}}[1..5];
	}
	my $aa3 = $aligns[2][$ai3];
	my $tname3 = align2::get_tname($aa3);
	my $qname3 = align2::get_qname($aa3);
	my $ovl1 = 0;
	my $ovl = 0;
	if ($tname3 eq $name1) {
		$ovl1 = czl_interval::czl_intersect_len( [@{$a1}[4..5]], [$tbeg3, $tend3] );
	}
	if ($ovl1>0) {
		$ovl = czl_interval::czl_intersect_len( [@{$a2}[4..5]], [$qbeg3, $qend3] );
	} else {
		if ($qname3 eq $name1) {
			$ovl1 = czl_interval::czl_intersect_len( [@{$a1}[4..5]], [$qbeg3, $qend3] );
			$ovl = czl_interval::czl_intersect_len( [@{$a2}[4..5]], [$tbeg3, $tend3] );
		}
	}
	if ($ovl>0) {
		$match_loci_count[2]++;
		$exon_stat[0]{$name}[$i1]{copy_type} = 2;
	} else {
		$match_loci_count[5]++;
		$exon_stat[0]{$name}[$i1]{copy_type} = 5;
	}
}
}

for (my $k=0; $k<8; $k++) {
	print STDERR "$k\t$match_loci_count[$k]\n";
}

$file = "${out_prefix}exon.triplet.txt";
$fout = base_io::czl_open($file, "w");

foreach my $name (keys(%{$exons[0]})) {
for (my $i1=0; $i1<@{$exons[0]{$name}}; $i1++) {
	my $stat='OK';
	my $type = $exon_stat[0]{$name}[$i1]{copy_type};
	if ($type>2) { next; }
	my $exon1 = $exons[0]{$name}[$i1];
	my ($a1,$aa1,$aid1,$name1,$id1);
	my ($a2,$aa2,$aid2,$name2,$id2);
	$id1 = ".";
	$id2 = ".";
	$a1 = $map_exons[0]{$name}[$i1][0];
	if (defined $a1) {
		$aa1 = $aligns[0][$a1->[0]];
		$aid1 = align2::get_id($aa1);
		$name1 = align2::get_tname($aa1);
		if (@{$a1->[8]}>0) {
			my @id1s;
			foreach my $ei (@{$a1->[8]}) {
				my $e = $exons[1]{$name1}[$ei];
				push @id1s, join(',', @$e[0..5]);
			}
			$id1 = join(';', @id1s);
		}
	}

	$a2 = $map_exons[1]{$name}[$i1][0];
	if (defined $a2) {
		$aa2 = $aligns[1][$a2->[0]];
		$aid2 = align2::get_id($aa2);
		$name2 = align2::get_tname($aa2);
		if (@{$a2->[8]}>0) {
			my @id2s;
			foreach my $ei (@{$a2->[8]}) {
				my $e = $exons[1]{$name2}[$ei];
				push @id2s, join(',', @$e[0..5]);
			}
			$id2 = join(';', @id2s);
		}
	}
	if (defined $a1) {
		if (defined $a2) {
		} else {
			$stat = 'LOST1';
		}
	} else {
		if (defined $a2) {
			$stat = 'LOST1';
		} else {
			$stat = 'LOST2';
		}
	}
	print $fout join("\t", @$exon1[0..5]), "\t";
	if (defined $a1) {
		print $fout join("\t", ($a1->[2],$a1->[3],$name1,$a1->[4], $a1->[5], $id1, 0, $a1->[6]), $aid1, $a1->[7]), "\t";
	} else {
		print $fout join("\t", ('.') x 10), "\t";
	}
	if (defined $a2) {
		print $fout join("\t", ($a2->[2],$a2->[3],$name2,$a2->[4], $a2->[5], $id2, 0, $a2->[6]), $aid2, $a2->[7]), "\t";
	} else {
		print $fout join("\t", ('.') x 10), "\t";
	}
	print $fout $stat, "\n";
}
}
close $fout;

$file = "${out_prefix}exon.sp2.txt";
$fout = base_io::czl_open($file, "w");

foreach my $name (keys(%{$exons[1]})) {
for (my $i1=0; $i1<@{$exons[1]{$name}}; $i1++) {
	my $exon1 = $exons[1]{$name}[$i1];
	my $a1s = $map_exons[2]{$name}[$i1];
	my $a2s = $map_exons[3]{$name}[$i1];
	my $a1;
	my $aa1;
	my $stat = 'OK';
	if (@$a1s==0 && @$a2s==0) {
		print $fout join("\t", @$exon1[0..5]), "\t";
		print $fout join("\t", (".") x 10), "\t";
		print $fout 'LOST1', "\n";
		next;
	} elsif (@$a1s>0 && @$a2s>0) {
		$stat ='Mul2';
		$a1 = $a1s->[0];
		$aa1 = $aligns[0][$a1->[0]];
	} elsif (@$a1s>1 || @$a2s>1) {
		$stat ='Mul1';
		if (@$a1s>1) {
			$a1 = $a1s->[0];
			$aa1 = $aligns[0][$a1->[0]];
		} else {
			$a1 = $a2s->[0];
			$aa1 = $aligns[1][$a1->[0]];
		}
	} elsif (@$a1s==1) {
		$a1 = $a1s->[0];
		$aa1 = $aligns[0][$a1->[0]];
	} elsif (@$a2s==1) {
		$a1 = $a2s->[0];
		$aa1 = $aligns[1][$a1->[0]];
	}
	my $aid1 = align2::get_id($aa1);
	my $name1 = align2::get_tname($aa1);
	my $name2 = align2::get_qname($aa1);
	my $id1 = ".";
	if (@{$a1->[8]}>0) {
		my @id1s;
		foreach my $ei (@{$a1->[8]}) {
			my $e = $exons[0]{$name2}[$ei];
			push @id1s, join(',', @$e[0..5]);
		}
		$id1 = join(';', @id1s);
	}
	print $fout join("\t", @$exon1[0..5]), "\t";
	print $fout join("\t", ($a1->[2],$a1->[3],$name2,$a1->[4], $a1->[5], $id1, 0, $a1->[6]), $aid1, $a1->[7]), "\t";
	print $fout $stat, "\n";
}
}
close $fout;


$file = "${out_prefix}exon.sp22.txt";
$fout = base_io::czl_open($file, "w");

foreach my $name (keys(%{$exons[1]})) {
for (my $i1=0; $i1<@{$exons[1]{$name}}; $i1++) {
	my $exon1 = $exons[1]{$name}[$i1];
	my $a1s = $map_exons[4]{$name}[$i1];
	my $a2s = $map_exons[5]{$name}[$i1];
	my $a1;
	my $aa1;
	my $stat='OK';
	my $name2;
	if (@$a1s==0 && @$a2s==0) {
		print $fout join("\t", @$exon1[0..5]), "\t";
		print $fout join("\t", (".") x 10), "\t";
		print $fout 'LOST2', "\n";
		next;
	} elsif (@$a1s>0 && @$a2s>0) {
		$a1 = $a1s->[0];
		$aa1 = $aligns[2][$a1->[0]];
		$name2 = align2::get_qname($aa1);
		$stat='Mul2';
	} elsif (@$a1s>1 || @$a2s>1) {
		if (@$a1s>0) {
			$a1 = $a1s->[0];
			$aa1 = $aligns[2][$a1->[0]];
			$name2 = align2::get_qname($aa1);
		} else {
			$a1 = $a2s->[0];
			$aa1 = $aligns[2][$a1->[0]];
			$name2 = align2::get_tname($aa1);
		}
		$stat='Mul1';
	} elsif (@$a1s==1) {
		$a1 = $a1s->[0];
		$aa1 = $aligns[2][$a1->[0]];
		$name2 = align2::get_qname($aa1);
	} elsif (@$a2s==1) {
		$a1 = $a2s->[0];
		$aa1 = $aligns[2][$a1->[0]];
		$name2 = align2::get_tname($aa1);
	}
	my $aid1 = align2::get_id($aa1);
	my $id1 = ".";
	if (@{$a1->[8]}>0) {
		my @id1s;
		foreach my $ei (@{$a1->[8]}) {
			my $e = $exons[1]{$name2}[$ei];
			push @id1s, join(',', @$e[0..5]);
		}
		$id1 = join(';', @id1s);
	}
	print $fout join("\t", @$exon1[0..5]), "\t";
	print $fout join("\t", ($a1->[2], $a1->[3], $name2, $a1->[4], $a1->[5], $id1, 0, $a1->[6]), $aid1, $a1->[7]);
	print $fout "\t", $stat, "\n";
}
}
close $fout;
