#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname fileparse);
use czl_io::base_io;
use czl_io::align2;
use czl_io::ucsc_chain;
use Set::IntervalTree;

sub usage()
{
print<<EOF;
Usage:
  $0 <-i1 GF_ZF.net.chain> <-i2 GF_GF.net.chain> <-o OUT>
Output:
EOF
}

my $fin;
my $fout;
my $file;
my $in_file;
my $in_file1;
my $in_file2;
my $out_file;
my $out_prefix;
my ($i,$j,$k);

if ($#ARGV<0) { usage(); exit 0; }

for (my $k=0; $k<=$#ARGV; $k++) {
	if ($ARGV[$k] =~ m/^(-i1)$/) {
		$in_file1 = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-i2)$/) {
		$in_file2 = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-o)$/) {
		$out_file = $ARGV[++$k];
	} elsif ($ARGV[$k] =~ m/^(-h|--help|-help)$/) {
		usage();
		exit 0;
	} else {
		die "No option '$ARGV[$k]'\n";
	}
}

my $aligns1 = ucsc_chain::load_all_align_from_chain($in_file1, 0);
@$aligns1 = sort { -(align2::get_score($a)<=>align2::get_score($b)) } @$aligns1;
my $aligns2;
$aligns2 = ucsc_chain::load_all_align_from_chain($in_file2, 0);
@$aligns2 = sort { -(align2::get_score($a)<=>align2::get_score($b)) } @$aligns2;

my @aligns = ($aligns1, $aligns2);
my @tree;
my @align_matches;

my ($tname, $tsize, $tbegin, $tend, $tstrand);
my ($qname, $qsize, $qbegin, $qend, $qstrand);

my $n1 = @{$aligns[0]};
my $n2 = @{$aligns[1]};
my @sel;
my %chr_part;
$sel[0] = [ (0) x $n1 ];
$sel[1] = [ (0) x $n2 ];
$k=0;
foreach my $k (0..1) {
	# k=0: species2 -- species1,  2:1
	# k=1: species2 -- species2,  1:1
for (my $ai=0; $ai<@{$aligns[$k]}; $ai++) {
	my $aa = $aligns[$k][$ai];
	($tname, $tsize, $tbegin, $tend, $tstrand) = align2::get_tcoord($aa);
	($qname, $qsize, $qbegin, $qend, $qstrand) = align2::get_qcoord($aa);
	my $regs = align2::czl_align2_to_interval_array2($aa, 3, 50);
	my $tregs = $regs->[0];
	my $qregs = $regs->[1];
	if (!exists $tree[2*$k+0]{$tname}) { $tree[2*$k+0]{$tname} = Set::IntervalTree->new; }
	if (!exists $tree[2*$k+1]{$qname}) { $tree[2*$k+1]{$qname} = Set::IntervalTree->new; }
	my $tree1=$tree[2*$k+0]{$tname};
	my $tree2=$tree[2*$k+1]{$qname};
#	if (!exists $tree[2]{$qname}) { $tree[2]{$qname} = Set::IntervalTree->new; }
#	if (!exists $tree[3]{$qname}) { $tree[3]{$qname} = Set::IntervalTree->new; }
	my $ql = czl_interval::czl_interval_len($qregs);
	my $tl = czl_interval::czl_interval_len($tregs);
	$align_matches[$k][$ai] = [$tl, $ql];
	my $ql_ovl = 0;
	my $tl_ovl = 0;
	if ($k==1 && $tname eq $qname) { next; }
#	if ($k==1 && $tname eq $qname) {
## check overlap of target and query
#		my $ovl = czl_interval::czl_intersect_interval_array_len($tregs, $qregs, 0);
#		if ($ovl>0) {
#			next;
#		}
#	}
# check Q overlap
	my @qregs_ovl; # like ZF 
	my @tregs_ovl; # like GF
	foreach my $reg (@$qregs) {
		my ($beg,$end) = @$reg[0..1];
		my $ai_v = $tree2->fetch($beg, $end); 
		if ($k==0 && @$ai_v>1) {
			my @regs1;
			foreach my $a (@$ai_v) {
				push @regs1, [ @$a[3..4] ];
			}
			my $pileup = czl_interval::czl_pileup_interval(\@regs1, 0, [$beg,$end]);
			my $l=0;
			foreach my $a (@$pileup) {
				if ($a->[2]>1) { $l+=$a->[1]-$a->[0]; }
			}
			if ($l>0) {
				$ql_ovl += $l;
#				push @qregs_ovl, $ai_v;
			}
		} elsif ($k==1 && @$ai_v>1) {
			$ql_ovl += $end-$beg;
		}
	}
	if ($ql_ovl>=$ql*0.1) {
		next;
	}

	foreach my $reg (@$tregs) {
		my ($beg,$end) = @$reg[0..1];
		my $ai_v = $tree1->fetch($beg, $end); 
		if (@$ai_v>0) {
			$tl_ovl += $end-$beg;
#			push @tregs_ovl, $ai_v;
		}
#		if (@$ai_v>1) {
#			my @regs1;
#			foreach my $a (@$ai_v) {
#				push @regs1, [ @$a[1..2] ];
#			}
#			my $pileup = czl_interval::czl_pileup_interval(\@regs1, 0, [$beg,$end]);
#			my $l=0;
#			foreach my $a (@$pileup) {
#				if ($a->[2]>1) { $l+=$a->[1]-$a->[0]; }
#			}
#			if ($l>0) {
#				$tl_ovl += $l;
#				push @tregs_ovl, $ai_v;
#			}
#		}
	}
	if ($tl_ovl>=$tl*0.1) {
		next;
	}
	
	for (my $j=0; $j<@$tregs; $j++) {
		my ($tbeg, $tend) = @{$tregs->[$j]}[0..1];
		my ($qbeg, $qend) = @{$qregs->[$j]}[0..1];
		$tree1->insert([$ai,$tbeg,$tend,$qbeg,$qend], $tbeg, $tend); 
		$tree2->insert([$ai,$tbeg,$tend,$qbeg,$qend], $qbeg, $qend); 
	}
	$sel[$k][$ai] = 1;
}
}

my %idxs;
for (my $ai=0; $ai<@{$aligns[$k]}; $ai++) {
	my $aa = $aligns[$k][$ai];
	my $qname = align2::get_qname($aa);
	if ($sel[0][$ai]) { push @{$idxs{$qname}}, $ai; }
}
my @tree2;
$k=0;
foreach my $qname (sort keys(%idxs)) {
	my $idxs1 = $idxs{$qname};
	@$idxs1 = sort {align2::get_qbegin($aligns[0][$a]) <=>align2::get_qbegin($aligns[0][$b]) } @$idxs1;
	$tree2[0]{$qname} = Set::IntervalTree->new;
	$tree2[1]{$qname} = Set::IntervalTree->new;
	my @tree3 = ($tree2[0]{$qname}, $tree2[1]{$qname});
	for (my $j=0; $j<@$idxs1; $j++) {
		my $ai = $idxs1->[$j];
		my $aa = $aligns[$k][$ai];
		($qname, $qsize, $qbegin, $qend, $qstrand) = align2::get_qcoord($aa);
		($tname, $tsize, $tbegin, $tend, $tstrand) = align2::get_tcoord($aa);
		my $regs = align2::czl_align2_to_interval_array2($aa, 3, 50);
		my $tregs = $regs->[0];
		my $qregs = $regs->[1];
		my $ql = czl_interval::czl_interval_len($qregs);
		my @ql_ovl = (0,0);
# check Q overlap
		my @qregs_ovl; # like ZF 
		for (my $k1=0; $k1<2; $k1++) {
			foreach my $reg (@$qregs) {
				my ($beg,$end) = @$reg[0..1];
				my $ai_v = $tree3[$k1]->fetch($beg, $end); 
				if (@$ai_v>0) {
					my @regs1;
					foreach my $a (@$ai_v) {
						push @regs1, [ @$a[3..4] ];
					}
					my $pileup = czl_interval::czl_pileup_interval(\@regs1, 0, [$beg,$end]);
					my $l=0;
					foreach my $a (@$pileup) {
						if ($a->[2]>0) { $l+=$a->[1]-$a->[0]; }
					}
					if ($l>0) {
						$ql_ovl[$k1] += $l;
					}
				}
			}
		}
		my $k1;
		$k1 = $ql_ovl[0]<=$ql_ovl[1] ? 0 : 1;
		for (my $j=0; $j<@$qregs; $j++) {
			my $qreg = $qregs->[$j];
			my $treg = $tregs->[$j];
			my ($qbeg,$qend) = @$qreg[0..1];
			my ($tbeg,$tend) = @$treg[0..1];
			$tree3[$k1]->insert([$ai,$tbeg,$tend,$qbeg,$qend], $qbeg, $qend); 
		}
		$sel[0][$ai] = $k1+1;
	}
}

my @count;
map { $count[0][$_]++; } @{$sel[0]};

my @filt=@{$sel[0]};
map {$filt[$_]= ($sel[0][$_]==1 ? 0 : 1); } 0..($n1-1);
$file = "${out_file}1a.chain";
ucsc_chain::store_all_align_to_chain($file, $aligns[0], \@filt);

map {$filt[$_]= ($sel[0][$_]==2 ? 0 : 1); } 0..($n1-1);
$file = "${out_file}1b.chain";
ucsc_chain::store_all_align_to_chain($file, $aligns[0], \@filt);

$file = "${out_file}1.filt.chain";
map {$filt[$_]= ($sel[0][$_]==0 ? 0 : 1); } 0..($n1-1);
ucsc_chain::store_all_align_to_chain($file, $aligns[0], \@filt);

map {$filt[$_]= ($sel[1][$_]!=0 ? 0 : 1); } 0..($n2-1);
$file = "${out_file}ab1.chain";
ucsc_chain::store_all_align_to_chain($file, $aligns[1], \@filt);

map {$filt[$_]= ($sel[1][$_]==0 ? 0 : 1); } 0..($n2-1);
$file = "${out_file}ab1.filt.chain";
ucsc_chain::store_all_align_to_chain($file, $aligns[1], \@filt);

#{
#	my $a = $tree[0]{'LG1'}->fetch(30128559,30128851);
#	my $b = $tree[0]{'LG26'}->fetch(8834758,8835091);
#
#	map { print join("\t", @$_), "\n"} @$a;
#	print "\n";
#	map { print join("\t", @$_), "\n"} @$b;
#}


my @tuples;
for (my $ai=0; $ai<@{$aligns[1]}; $ai++) {
	if ($sel[1][$ai]==0) { next; }
	my $aa = $aligns[1][$ai];
	my ($name1, $size1, $begin1, $end1, $strand1) = align2::get_tcoord($aa);
	my ($name2, $size2, $begin2, $end2, $strand2) = align2::get_qcoord($aa);
	my $span1 = $end1-$begin1;
	my $span2 = $end2-$begin2;
	my $regs = align2::czl_align2_to_interval_array2($aa, 3, 50);
	my $regs1 = $regs->[0];
	my $regs2 = $regs->[1];
	my $l1 = czl_interval::czl_interval_len($regs1);
	my $l2 = czl_interval::czl_interval_len($regs2);
	my $l3 = 0;
	my $l_ovl1 = 0;
	my $l_ovl2 = 0;
	my @l_ovl = (0,0);
# check Q overlap
	my @regs_ovl1; # like ZF 
	my @regs_ovl2; # like GF
	my @l_ovl_by_ai;
	my @begs3=([],[]);
	my @ends3=([],[]);
	my @ais;
	for (my $j=0; $j<@$regs1; $j++) {
		my @begs;
		my @ends;
		$begs[0] = $regs1->[$j][0];
		$ends[0] = $regs1->[$j][1];
		$begs[1] = $regs2->[$j][0];
		$ends[1] = $regs2->[$j][1];
		my @ai_v;
		my @v;
		foreach my $k1 (0..1) {
			my $name = $k1==0 ? $name1 : $name2;
			if (!exists $tree[0]{$name}) {$ai_v[$k1] = []; }
			else { $ai_v[$k1] = $tree[0]{$name}->fetch($begs[$k1], $ends[$k1]);}
			map {
				my $ai3 = $_->[0];
				my $k3 = $sel[0][$ai3];
				if (exists $v[$k3-1][$k1]{$ai3}) {
					if ($_->[3] < $v[$k3-1][$k1]{$ai3}[0]) {
						$v[$k3-1][$k1]{$ai3}[0] = $_->[3];
					}
					if ($_->[4] > $v[$k3-1][$k1]{$ai3}[1]) {
						$v[$k3-1][$k1]{$ai3}[1] = $_->[4];
					}
				} else {
					$v[$k3-1][$k1]{$ai3}[0] = $_->[3];
					$v[$k3-1][$k1]{$ai3}[1] = $_->[4];
				}
			} @{$ai_v[$k1]};
		}
		foreach my $k1 (0..1) {
		foreach my $k2 (0..1) {
			my $v1 = $v[$k1][$k2];
			if (!defined $v1) { $v[$k1][$k2]=[]; next;}
			my @v2;
			foreach my $ai (keys(%$v1)) {
				push @v2, [@{$v1->{$ai}}, $ai];
			}
			$v[$k1][$k2] = \@v2;
		}
		}
#		foreach my $k1 (0..1) {
#		foreach my $k2 (0..1) {
#			my $v1 = $v[$k1][$k2];
#			if (@$v1==0) { next; }
#			my ($c,$to_c) = czl_interval::czl_group_by_dist($v1, 0);
#			my $u = czl_interval::czl_merge_interval_by_group($v1, $c);
#			$v[$k1][$k2] = $u;
#		}
#		}
		foreach my $k1 (0..1) {
			my $u1 = $v[$k1][0];
			my $u2 = $v[1-$k1][1];
			my $sorted=0;
			my $u3 = czl_interval::czl_intersect_interval_array($u1,$u2,$sorted);
			if (@$u3>0) {
				my $l = czl_interval::czl_interval_len($u3);
				foreach my $k2 (0..1) {
					if (!defined $begs3[$k1][$k2] || $begs[$k2]<$begs3[$k1][$k2]) {
						$begs3[$k1][$k2] = $begs[$k2];
					}
					if (!defined $ends3[$k1][$k2] || $ends[$k2]>$ends3[$k1][$k2]) {
						$ends3[$k1][$k2] = $ends[$k2];
					}
				}
				$l_ovl[$k1] += $l;

				foreach my $reg (@$u3) {
					my $j1 = $reg->[2];
					my $j2 = $reg->[3];
					my $ai1 = $u1->[$j1][2];
					my $ai2 = $u2->[$j2][2];
					$ais[$k1]{$ai1}{$ai2}+=$reg->[1]-$reg->[0];
				}
			}
		}
#		foreach my $k1 (0..1) {
#			my $u1 = $v[$k1][0];
#			my $u2 = $v[1-$k1][1];
#			foreach my $ai1 (keys(%$u1)) {
#				my $u3 = $u1->{$ai1};
#				foreach my $ai2 (keys(%$u2)) {
#					my $u4 = $u2->{$ai2};
#					my $sorted=0;
#					my $l = czl_interval::czl_intersect_interval_array_len($u3,$u4,$sorted);
#					if ($l>0) { $l_ovl_by_ai[$k1]{$ai1}{$ai2}+=$l; }
#				}
#			}
#		}
	}
#	foreach my $k1 (0..1) {
#		foreach my $ai1 (keys(%{$l_ovl_by_ai[$k1]})) {
#			foreach my $ai2 (keys(%{$l_ovl_by_ai[$k1]{$ai1}})) {
#				my $l = $l_ovl_by_ai[$k1]{$ai1};
#				my $l1 = $align_matches[0][$ai1];
#				my $l2 = $align_matches[0][$ai2];
#			}
#		}
#	}
	
	my @span3;
	foreach my $k1 (0..1) {
	foreach my $k2 (0..1) {
		if (!defined $begs3[$k1][$k2]) { $span3[$k1][$k2]=0; }
		else {$span3[$k1][$k2] = $ends3[$k1][$k2]-$begs3[$k1][$k2]; }
	}
	}
	my $k1 = $span3[0][0]+$span3[0][1] >= $span3[1][0]+$span3[1][1] ? 0 : 1;
	if ($span3[$k1][0]<$span1*0.25 && $span3[$k1][1]<$span2*0.25) {
		$sel[1][$ai] = 0;
		next;
	}

	$sel[1][$ai] = $k1+1;

	if ($k1==1) {
		align2::swap_QT($aa);
		if (align2::get_tstrand($aa) eq '-') {
			align2::flip($aa);
		}
	}

	foreach my $ai1 (keys(%{$ais[$k1]})) {
	foreach my $ai2 (keys(%{$ais[$k1]{$ai1}})) {
		push @tuples, [$ai, $ai1, $ai2, $ais[$k1]{$ai1}{$ai2}];
	}
	}
}
map { $count[1][$_]++; } @{$sel[1]};


map {$filt[$_]= ($sel[1][$_]!=0 ? 0 : 1); } 0..($n2-1);
$file = "${out_file}ab2.chain";
ucsc_chain::store_all_align_to_chain($file, $aligns[1], \@filt);

map {$filt[$_]= ($sel[1][$_]==0 ? 0 : 1); } 0..($n2-1);
$file = "${out_file}ab2.filt.chain";
ucsc_chain::store_all_align_to_chain($file, $aligns[1], \@filt);


$file = "${out_file}chain_triplet.txt";
$fout = base_io::czl_open($file, "w");
foreach my $a (@tuples) {
	print $fout align2::get_id($aligns[0][$a->[1]]), "\t";
	print $fout align2::get_id($aligns[0][$a->[2]]), "\t";
	print $fout align2::get_id($aligns[1][$a->[0]]), "\t";
	print $fout $a->[3], "\n";
}
close $fout;
