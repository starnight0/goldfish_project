#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Cwd qw(abs_path getcwd);
use File::Basename qw(basename dirname fileparse);
use czl_common;
use czl_io::base_io;
use czl_io::align2;
use czl_math;
#use Set::IntervalTree;

sub usage()
{
print<<EOF;
Usage:
  $0 <-i INFILE> <-o OUT> [--iden_thres1 IDEN1] [--iden_thres2 IDEN2]
Output:
EOF
}

my ($i,$j,$k);
my $fin;
my $fout;
my $file;
my $in_file;
my $out_file;
my $out_prefix;
my %conf = (
	in_file  => [ ['in_file'] , 'i', 1, 1, [0], undef],
	out_file => [ ['out_file'], 'o', 1, 1, [0], undef],
	overlap_frac_thres => [ ['overlap_frac_thres'], 'ovlf', 1, 0, [0], 0.1],
	iden_thres1 => [ [], undef, 1, 0, [0], 0.8],
	iden_thres2 => [ [], undef, 1, 0, [0], 0.75],
		);

if ($#ARGV<0) { usage(); exit 0; }

czl_common::czl_parse_conf(\@ARGV, \%conf);
my $iden_thres1 = $conf{iden_thres1}[czl_common::CONF_VALUE];
my $iden_thres2 = $conf{iden_thres2}[czl_common::CONF_VALUE];
my $tcopy = 1;
my $qcopy = 2;
my $ovlf_thres = $conf{overlap_frac_thres}[czl_common::CONF_VALUE];
my $ovlf_high_thres = 0.9;
$out_file = $conf{out_file}[czl_common::CONF_VALUE];


my $align_v = align2::load_from_psl($conf{in_file}[czl_common::CONF_VALUE], 0);

if (@$align_v==0) { die("No alignments\n");  }

# @$align_v = sort { -(align2::get_iden_n($a)<=>align2::get_iden_n($b)) } @$align_v;

my @align;
my $n=0;
my $i0=0;
for (my $i=0; $i<@$align_v; $i++) {
	my $aa = $align_v->[$i];
#    my $tname = align2::get_tname($aa);
#    my $qname = align2::get_qname($aa);
#    my $tbegin = align2::get_tbegin($aa);
#    my $tend = align2::get_tend($aa);
    my $iden_n = align2::get_iden_n($aa);
    align2::set_score($aa, $iden_n);
    my $mis_match = align2::get_mismatch($aa);
#    my $t_ins_bp = align2::get_unaligned_l($aa,0);
#    my $q_ins_bp = align2::get_unaligned_l($aa,1);
#    my $alen  = align2::get_align_len($aa);
#    my $ins_bp = $t_ins_bp < $q_ins_bp ? $t_ins_bp : $q_ins_bp ;
    my $iden1 = $iden_n/($iden_n+$mis_match);
#    my $iden2 = $iden_n/($iden_n+$mis_match+$ins_bp);
    if ($iden1>=$iden_thres1) {
		if ($i0!=$i) { $align_v->[$i0] = $align_v->[$i]; }
		$i0++;
	}
#	push @{$align[0]{$tname}{$qname}}, $aa;
#	push @{$align[1]{$qname}{$tname}}, $aa;
#    if ($iden2<$iden_thres2) { next; }
}
$#$align_v = $i0;
my $na = @$align_v;

#my $ttree = Set::IntervalTree->new;
#my $qtree = Set::IntervalTree->new;

my @align_is1;
my @edges;
my %edge;
my %idx;
for (my $i=0; $i<@$align_v; $i++) {
	my $aa = $align_v->[$i];
	my $tname  = align2::get_tname($aa);
	my $tbegin = align2::get_tbegin($aa);
	push @{$idx{$tname}}, [$tbegin, $i];
}
my $m=0;
my @high_ovl;
foreach my $tname (keys(%idx)) {
	my $idx1 = $idx{$tname};
	@$idx1 = sort {$a->[0] <=> $b->[0]} @$idx1;
	if (@$idx1==1) { next; }
for (my $j0=0; $j0<@$idx1; $j0++) {
	my $i0 = $idx1->[$j0][1];
	my $aa0 = $align_v->[$i0];
	my $beg0 = align2::get_tbegin($aa0);
	my $end0 = align2::get_tend($aa0);
	my $l0 = $end0-$beg0;
	for (my $j=$j0+1; $j<@$idx1; $j++) {
		my $i = $idx1->[$j][1];
		my $aa = $align_v->[$i];
		my $beg = align2::get_tbegin($aa);
		my $end = align2::get_tend($aa);
		my $l = $end-$beg;
		my ($min_l, $max_l);
		if ($l0<$l) {$min_l = $l0; $max_l=$l; }
		else {$min_l = $l; $max_l=$l0; }
		if ($beg<$end0) {
			$m++;
			my $ovl;
			if ($end<$end0) {$ovl = $l;}
			else {$ovl = $end0-$beg;}
			if ($ovl >= $max_l * $ovlf_high_thres) {
				$high_ovl[$i0]++;
				if ($high_ovl[$i0]>20) { last; }
			}
			if ($ovl >= $min_l * $ovlf_thres) {
				push @edges, [$i0, $i];
			}
		} else {
			last;
		}
	}
}
}

undef %idx;
for (my $i=0; $i<@$align_v; $i++) {
	if (defined $high_ovl[$i]) { next; }
	my $aa = $align_v->[$i];
	my $qname  = align2::get_qname($aa);
	my $qbegin = align2::get_qbegin($aa);
	push @{$idx{$qname}}, [$qbegin, $i];
}
foreach my $qname (keys(%idx)) {
	my $idx1 = $idx{$qname};
	@$idx1 = sort {$a->[0] <=> $b->[0]} @$idx1;
	if (@$idx1==1) { next; }
for (my $j0=0; $j0<@$idx1; $j0++) {
	my $i0 = $idx1->[$j0][1];
	my $aa0 = $align_v->[$i0];
	my $beg0 = align2::get_qbegin($aa0);
	my $end0 = align2::get_qend($aa0);
	my $l0 = $end0-$beg0;
	for (my $j=$j0+1; $j<@$idx1; $j++) {
		my $i = $idx1->[$j][1];
		my $aa = $align_v->[$i];
		my $beg = align2::get_qbegin($aa);
		my $end = align2::get_qend($aa);
		my $l = $end-$beg;
		my $min_l = $l0<$l ? $l0 : $l;
		my $max_l = $l0>$l ? $l0 : $l;
		if ($beg<$end0) {
			my $ovl;
			if ($end<$end0) {$ovl = $l;}
			else {$ovl = $end0-$beg;}
			if ($ovl >= $max_l * $ovlf_high_thres) {
				$high_ovl[$i0]++;
				if ($high_ovl[$i0]>20) { last; }
			}
			if ($ovl >= $min_l * $ovlf_thres) {
				push @edges, [$i0, $i];
			}
		} else {
			last;
		}
	}
}
}
undef %idx;
@edges = sort {$a->[0]<=>$b->[0] or $a->[1]<=>$b->[1] } @edges;
$i = 0;
for (my $j=0; $j<@edges; $j++) {
	my $f = 0;
	for (my $k=0; $k<2; $k++) {
		if (defined $high_ovl[$edges[$j][$k]]) { $f=1; last; }
	}
	if ($i>0 && $edges[$j][0]==$edges[$i-1][0] && $edges[$j][1]==$edges[$i-1][1]) { $f=1; }
	if (!$f) {
		if ($i!=$j) { $edges[$i]=$edges[$j];}
		$i++;
	}
}
$#edges = $i-1;


my ($c, $to_c, $ungroup) = czl_math::single_link_cluster1(\@edges);
my $is_group;
my @selected;
my @c_edge;
my @align_idxs;
my @is_sel;
foreach my $e (@edges) {
	my $ci = $to_c->{$e->[0]};
	push @{$c_edge[$ci]}, $e;
}
foreach my $cid (0..$#$c) {
	my $c1 = $c->[$cid];
	@$c1=sort {$a<=>$b} @$c1;
	my %vw;
	map {$vw{$_} = align2::get_iden_n($align_v->[$_])} @$c1;
#	my ($set, $w) = czl_math::max_weight_independent_set_greed(\%vw, $c_edge[$cid]);
	my ($set, $w) = czl_math::max_weight_independent_set_local_search(\%vw, $c_edge[$cid], 20, 1000, 20, 29384, 1);
	if (!czl_math::is_independent_set($c_edge[$cid], $set)) { die; } 
	push @{$align_idxs[0]}, @$set;
	foreach my $i (@$set) { $is_sel[$i]=0; }
	
}

valid1($align_idxs[0]);

my @trank;
my @qrank;
my @idx_by_tbegin;
my @idx_by_qbegin;
for (my $i=0; $i<@$align_v; $i++) {
	my $aa = $align_v->[$i];
	my $tname = align2::get_tname($aa);
	my $tbegin = align2::get_tbegin($aa);
	push @idx_by_tbegin, [$i, $tname, $tbegin];
	my $qname = align2::get_qname($aa);
	my $qbegin = align2::get_qbegin($aa);
	push @idx_by_qbegin, [$i, $qname, $qbegin];
}
@idx_by_tbegin = sort {$a->[1] cmp $b->[1] or $a->[2]<=>$b->[2]} @idx_by_tbegin;
@idx_by_qbegin = sort {$a->[1] cmp $b->[1] or $a->[2]<=>$b->[2]} @idx_by_qbegin;
#my @rank_by_tbegin = (0) x $na;
#my @rank_by_qbegin = (0) x $na;
#for (my $j=0; $j<=@idx_by_tbegin; $j++) {
#	my $i = $idx_by_tbegin[$j][0];
#	$rank_by_tbegin[$i] = $j;
#}
#for (my $j=0; $j<=@idx_by_qbegin; $j++) {
#	my $i = $idx_by_qbegin[$j][0];
#	$rank_by_qbegin[$i] = $j;
#}

my @filt=(0) x $na;
my $nf = 0;
for (my $j0=0; $j0<@idx_by_tbegin; $j0++) {
	my $i0 = $idx_by_tbegin[$j0][0];
	if ($filt[$i0]==1) { next; }
	my $aa0 = $align_v->[$i0];
	my $name0 = align2::get_tname($aa0);
	my $beg0 = align2::get_tbegin($aa0);
	my $end0 = align2::get_tend($aa0);
	my $l0 = $end0-$beg0;
	if (defined $is_sel[$i0]) {
		for (my $j=$j0+1; $j<@idx_by_tbegin; $j++) {
			my $i = $idx_by_tbegin[$j][0];
			if ($filt[$i]==1) { next; }
			my $aa = $align_v->[$i];
			my $name = align2::get_tname($aa);
			my $beg = align2::get_tbegin($aa);
			my $end = align2::get_tend($aa);
			my $l = $end-$beg;
			my ($min_l, $max_l);
			if ($l0<$l) {$min_l = $l0; $max_l=$l; }
			else {$min_l = $l; $max_l=$l0; }
			my $is_ovl = 0;
			if ($name eq $name0 && $beg<$end0) {
				$m++;
				my $ovl;
				if ($end<$end0) {$ovl = $l;}
				else {$ovl = $end0-$beg;}
				if ($ovl >= $min_l * $ovlf_thres) { $is_ovl = 1; }
			}
			if ($is_ovl) {
				if (defined $is_sel[$i]) {
					die;
				}
				$filt[$i]=1;
				$nf++;
			} else {
				$j0=$j-1;
				last;
			}
		}
	}
}
for (my $j0=0; $j0<@idx_by_qbegin; $j0++) {
	my $i0 = $idx_by_qbegin[$j0][0];
	if ($filt[$i0]==1) { next; }
	my $aa0 = $align_v->[$i0];
	my $name0 = align2::get_qname($aa0);
	my $beg0 = align2::get_qbegin($aa0);
	my $end0 = align2::get_qend($aa0);
	my $l0 = $end0-$beg0;
	if (defined $is_sel[$i0]) {
		for (my $j=$j0+1; $j<@idx_by_qbegin; $j++) {
			my $i = $idx_by_qbegin[$j][0];
			if ($filt[$i]==1) { next; }
			my $aa = $align_v->[$i];
			my $name = align2::get_qname($aa);
			my $beg = align2::get_qbegin($aa);
			my $end = align2::get_qend($aa);
			my $l = $end-$beg;
			my ($min_l, $max_l);
			if ($l0<$l) {$min_l = $l0; $max_l=$l; }
			else {$min_l = $l; $max_l=$l0; }
			my $is_ovl = 0;
			if ($name eq $name0 && $beg<$end0) {
				$m++;
				my $ovl;
				if ($end<$end0) {$ovl = $l;}
				else {$ovl = $end0-$beg;}
				if ($ovl >= $min_l * $ovlf_thres) { $is_ovl = 1; }
			}
			if ($is_ovl) {
				if (defined $is_sel[$i]) {
					die;
				}
				$filt[$i]=1;
				$nf++;
			} else {
				$j0=$j-1;
				last;
			}
		}
	}
}
@edges = ();
my ($j0, $j1);
for ($j0=0; $j0<@idx_by_tbegin; $j0++) {
	my $i0 = $idx_by_tbegin[$j0][0];
	if ($filt[$i0]==1) { next; }
	if (defined $is_sel[$i0]) { next; }
	my $aa0 = $align_v->[$i0];
	my $name0 = align2::get_tname($aa0);
	my $beg0 = align2::get_tbegin($aa0);
	my $end0 = align2::get_tend($aa0);
	my $l0 = $end0-$beg0;
	if (!defined $j1 || $j1<=$j0) { $j1=$j0+1; }
	my $j;
	for ($j=$j1; $j<@idx_by_tbegin; $j++) {
		my $i = $idx_by_tbegin[$j][0];
		if ($filt[$i]==1) { next; }
		my $aa = $align_v->[$i];
		my $name = align2::get_tname($aa);
		my $beg = align2::get_tbegin($aa);
		my $end = align2::get_tend($aa);
		my $l = $end-$beg;
		my ($min_l, $max_l);
		if ($l0<$l) {$min_l = $l0; $max_l=$l; }
		else {$min_l = $l; $max_l=$l0; }

		my $is_ovl = 0;
		if ($name eq $name0 && $beg<$end0) {
			$m++;
			my $ovl;
			if ($end<$end0) {$ovl = $l;}
			else {$ovl = $end0-$beg;}
			if ($ovl >= $min_l * $ovlf_thres) { $is_ovl = 1; }
		}
		if ($is_ovl) {
			if (defined $is_sel[$i]) { $filt[$i0]=1; last; }
			push @edges, [$i0,$i];
		} else {
			$j1 = $j;
			last;
		}
	}
	if ($j1==$j0+1) { undef $j1; }
}
undef $j1;
for ($j0=0; $j0<@idx_by_qbegin; $j0++) {
	my $i0 = $idx_by_qbegin[$j0][0];
	if ($filt[$i0]==1) { next; }
	if (defined $is_sel[$i0]) { next; }
	my $aa0 = $align_v->[$i0];
	my $name0 = align2::get_qname($aa0);
	my $beg0 = align2::get_qbegin($aa0);
	my $end0 = align2::get_qend($aa0);
	my $l0 = $end0-$beg0;
	if (!defined $j1 || $j1<=$j0) { $j1=$j0+1; }
	my $j;
	for ($j=$j1; $j<@idx_by_qbegin; $j++) {
		my $i = $idx_by_qbegin[$j][0];
		if ($filt[$i]==1) { next; }
		my $aa = $align_v->[$i];
		my $name = align2::get_qname($aa);
		my $beg = align2::get_qbegin($aa);
		my $end = align2::get_qend($aa);
		my $l = $end-$beg;
		my ($min_l, $max_l);
		if ($l0<$l) {$min_l = $l0; $max_l=$l; }
		else {$min_l = $l; $max_l=$l0; }

		my $is_ovl = 0;
		if ($name eq $name0 && $beg<$end0) {
			$m++;
			my $ovl;
			if ($end<$end0) {$ovl = $l;}
			else {$ovl = $end0-$beg;}
			if ($ovl >= $min_l * $ovlf_thres) { $is_ovl = 1; }
		}
		if ($is_ovl) {
			if (defined $is_sel[$i]) { $filt[$i0]=1; last; }
			push @edges, [$i0,$i];
		} else {
			$j1 = $j;
			last;
		}
	}
	if ($j1==$j0+1) { undef $j1; }
}
@edges = sort {$a->[0]<=>$b->[0] or $a->[1]<=>$b->[1] } @edges;
$i = 0;
for (my $j=0; $j<@edges; $j++) {
	my $e = $edges[$j];
	my $f = 0;
	if ($filt[$e->[0]]==1 || $filt[$e->[1]]==1) { $f=1; }
	if ($i>0 && $edges[$j][0]==$edges[$i-1][0] && $edges[$j][1]==$edges[$i-1][1]) { $f=1; }
	if (!$f) {
		if ($i!=$j) { $edges[$i]=$edges[$j];}
		$i++;
		$filt[$e->[0]]=2;
		$filt[$e->[1]]=2;
	}
}
$#edges = $i-1;

($c, $to_c, $ungroup) = czl_math::single_link_cluster1(\@edges);

#@c_edge=();
#foreach my $e (@edges) {
#	my $ci = $to_c->{$e->[0]};
#	push @{$c_edge[$ci]}, $e;
#}
#foreach my $cid (0..$#$c) {
#	my $c1 = $c->[$cid];
#	@$c1=sort {$a<=>$b} @$c1;
#	my %vw;
#	map {$vw{$_} = align2::get_iden_n($align_v->[$_])} @$c1;
#	my ($set, $w) = czl_math::max_weight_independent_set_greed(\%vw, $c_edge[$cid]);
#	my ($set, $w) = czl_math::max_weight_independent_set_local_search(\%vw, $c_edge[$cid], 20, 1000, 20, 29384, 1);
#	if (!czl_math::is_independent_set($c_edge[$cid], $set)) { die; } 
#	push @{$align_idxs[0]}, @$set;
#	foreach my $i (@$set) { $is_sel[$i]=0; }
#	valid1($align_idxs[0]);
#}


for (my $cid=0; $cid<@$c; $cid++) {
	my $c1 = $c->[$cid];
	@$c1 = sort { -(align2::get_score($align_v->[$a]) <=> align2::get_score($align_v->[$b])) } @$c1;
	my ($j0, $j1);
	my @set = ($c1->[0]);
	for ($j1=0; $j1<@$c1; $j1++) {
		my $i1 = $c1->[$j1];
		my $aa1 = $align_v->[$i1];
		my $name1 = align2::get_tname($aa1);
		my $beg1 = align2::get_tbegin($aa1);
		my $end1 = align2::get_tend($aa1);
		my $l1 = $end1-$beg1;
		my $is_ovl=0;
		for ($j0=0; $j0<@set; $j0++) {
			my $i0 = $set[$j0];
			my $aa0 = $align_v->[$i0];
			my $name0 = align2::get_tname($aa0);
			my $beg0 = align2::get_tbegin($aa0);
			my $end0 = align2::get_tend($aa0);
			my $l0 = $end0-$beg0;
			my $min_l = $l0<$l1 ? $l0 : $l1;
			my $b = $beg0 > $beg1 ? $beg0 : $beg1;
			my $e = $end0 < $end1 ? $end0 : $end1;
			my $ovl = $e-$b;
			if ($ovl >= $min_l * $ovlf_thres) {
				$is_ovl = 1; last;
			}
		}
		if ($is_ovl) { next; }
		$name1 = align2::get_qname($aa1);
		$beg1 = align2::get_qbegin($aa1);
		$end1 = align2::get_qend($aa1);
		for ($j0=0; $j0<@set; $j0++) {
			my $i0 = $set[$j0];
			my $aa0 = $align_v->[$i0];
			my $name0 = align2::get_qname($aa0);
			my $beg0 = align2::get_qbegin($aa0);
			my $end0 = align2::get_qend($aa0);
			my $l0 = $end0-$beg0;
			my $min_l = $l0<$l1 ? $l0 : $l1;
			my $b = $beg0 > $beg1 ? $beg0 : $beg1;
			my $e = $end0 < $end1 ? $end0 : $end1;
			my $ovl = $e-$b;
			if ($ovl >= $min_l * $ovlf_thres) {
				$is_ovl = 1; last;
			}
		}
		if ($is_ovl) { next; }
		push @set, $i1;
	}
	push @{$align_idxs[0]}, @set;
	map { $is_sel[$_]=0; } @set;
	valid1($align_idxs[0]);
}

for (my $i=0; $i<$na; $i++) {
	if (defined $is_sel[$i]) { next; }
	if ($filt[$i]) { next; }
	push @{$align_idxs[0]}, $i;
	$is_sel[$i]=0;
	valid1($align_idxs[0]);
}

foreach my $i (sort { (align2::get_tname($align_v->[$a]) cmp align2::get_tname($align_v->[$b])) or align2::get_tbegin($align_v->[$a]) <=> align2::get_tbegin($align_v->[$b])} @{$align_idxs[0]}) {
	my $aa = $align_v->[$i];
	my $tname = align2::get_tname($aa);
	my $tbegin = align2::get_tbegin($aa);
	my $tend = align2::get_tend($aa);
	my $qname = align2::get_qname($aa);
	my $qbegin = align2::get_qbegin($aa);
	my $qend = align2::get_qend($aa);
	print join("\t", ($tname,$tbegin,$tend, $qname,$qbegin,$qend)), "\n";
}


valid1($align_idxs[0]);
valid2($align_idxs[0]);

sub valid1
{
	my $idxs = shift;
	my @idxs1 = sort {align2::get_tname($align_v->[$a]) cmp align2::get_tname($align_v->[$b]) or align2::get_tbegin($align_v->[$a]) <=> align2::get_tbegin($align_v->[$b]) } @$idxs;
	for (my $j0=0; $j0<$#idxs1; $j0++) {
		my $i0 = $idxs1[$j0];
		my $aa0 = $align_v->[$i0];
		my $name0 = align2::get_tname($aa0);
		my $beg0 = align2::get_tbegin($aa0);
		my $end0 = align2::get_tend($aa0);
		my $l0 = $end0-$beg0;
		my $j1 = $j0+1;
		my $i1 = $idxs1[$j1];
		my $aa1 = $align_v->[$i1];
		my $name1 = align2::get_tname($aa1);
		if ($name0 ne $name1) { next; }
		my $beg1 = align2::get_tbegin($aa1);
		my $end1 = align2::get_tend($aa1);
		my $l1 = $end1-$beg1;
		my $min_l = $l0<$l1 ? $l0 : $l1;
		my $b = $beg0 > $beg1 ? $beg0 : $beg1;
		my $e = $end0 < $end1 ? $end0 : $end1;
		my $ovl = $e-$b;
		if ($ovl >= $min_l * $ovlf_thres) {
			die;
		}
	}
	@idxs1 = sort {align2::get_qname($align_v->[$a]) cmp align2::get_qname($align_v->[$b]) or align2::get_qbegin($align_v->[$a]) <=> align2::get_qbegin($align_v->[$b]) } @$idxs;
	for (my $j0=0; $j0<$#idxs1; $j0++) {
		my $i0 = $idxs1[$j0];
		my $aa0 = $align_v->[$i0];
		my $name0 = align2::get_qname($aa0);
		my $beg0 = align2::get_qbegin($aa0);
		my $end0 = align2::get_qend($aa0);
		my $l0 = $end0-$beg0;
		my $j1 = $j0+1;
		my $i1 = $idxs1[$j1];
		my $aa1 = $align_v->[$i1];
		my $name1 = align2::get_qname($aa1);
		if ($name0 ne $name1) { next; }
		my $beg1 = align2::get_qbegin($aa1);
		my $end1 = align2::get_qend($aa1);
		my $l1 = $end1-$beg1;
		my $min_l = $l0<$l1 ? $l0 : $l1;
		my $b = $beg0 > $beg1 ? $beg0 : $beg1;
		my $e = $end0 < $end1 ? $end0 : $end1;
		my $ovl = $e-$b;
		if ($ovl >= $min_l * $ovlf_thres) {
			die;
		}
	}
}

sub valid2
{
	my $idxs = shift;
	my $na = $#$align_v+1;
	my @is_sel = (0) x $na;
	map { $is_sel[$_] = 0; } @$idxs;
	for (my $i0=0; $i0<$na; $i0++) {
		if (defined $is_sel[$i0]) { next; }
		my $is_ovl=0;
		my $aa0 = $align_v->[$i0];
		my $name0 = align2::get_tname($aa0);
		my $beg0 = align2::get_tbegin($aa0);
		my $end0 = align2::get_tend($aa0);
		my $l0 = $end0-$beg0;
		for (my $j1=0; $j1<$#$idxs; $j1++) {
			my $i1 = $idxs->[$j1];
			my $aa1 = $align_v->[$i1];
			my $name1 = align2::get_tname($aa1);
			if ($name0 ne $name1) { next; }
			my $beg1 = align2::get_tbegin($aa1);
			my $end1 = align2::get_tend($aa1);
			my $l1 = $end1-$beg1;
			my $min_l = $l0<$l1 ? $l0 : $l1;
			my $b = $beg0 > $beg1 ? $beg0 : $beg1;
			my $e = $end0 < $end1 ? $end0 : $end1;
			my $ovl = $e-$b;
			if ($ovl >= $min_l * $ovlf_thres) {
				$is_ovl = 1;
			}
		}
		$name0 = align2::get_qname($aa0);
		$beg0 = align2::get_qbegin($aa0);
		$end0 = align2::get_qend($aa0);
		$l0 = $end0-$beg0;
		for (my $j1=0; $j1<$#$idxs; $j1++) {
			my $i1 = $idxs->[$j1];
			my $aa1 = $align_v->[$i1];
			my $name1 = align2::get_qname($aa1);
			if ($name0 ne $name1) { next; }
			my $beg1 = align2::get_qbegin($aa1);
			my $end1 = align2::get_qend($aa1);
			my $l1 = $end1-$beg1;
			my $min_l = $l0<$l1 ? $l0 : $l1;
			my $b = $beg0 > $beg1 ? $beg0 : $beg1;
			my $e = $end0 < $end1 ? $end0 : $end1;
			my $ovl = $e-$b;
			if ($ovl >= $min_l * $ovlf_thres) {
				$is_ovl = 1;
			}
		}
		if (!$is_ovl) {
			die;
		}
	}
}
