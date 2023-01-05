use strict;
use warnings;
use Data::Dumper;

my ($i,$j,$k);
my $n = 0;
my $copy = 2;
my @gene; # sp1, sp2
my @exon; # sp1, sp2
my $triple_file = $ARGV[0];
my $sp2_file    = $ARGV[1];
my $sp22_file   = $ARGV[2];
my $out_prefix = $ARGV[3];
open IN, "<$triple_file" or die;
while(<IN>) {
	if (m/^\s*$/) { next; }
	s/\s+$//;
	my @t = split /\t/;
	if ($t[26]=~m/Mul/i) { next; }
	my @ids;
	my @begs;
	my @ends;
	$begs[0] = $t[1];
	$ends[0] = $t[2];
	$ids[0] = $t[3];
	$begs[1] = $t[9];
	$ends[1] = $t[10];
	$ids[1] = $t[11];
	$begs[2] = $t[19];
	$ends[2] = $t[20];
	$ids[2] = $t[21];
	my @idss;
	my @gids;
	my @ee = split /:/, $ids[0];
	my $gid0='.';
	my $is_exon; # if mapped to an exon, set to 1;
	if ($ids[0] =~ m/exon/i) { $is_exon=1; $gid0=$ee[5];}
	else { $is_exon=0; }
	my $gid1 = '.';
	my $l1=0;
	my $gid2 = '.';
	my $l2=0;
	my @lens = (0,0,0);
	my @exon1 = @t[0..5];
	$exon1[8] = $t[26];
	$lens[0] = $t[2]-$t[1];
	for (my $k=1; $k<1+$copy; $k++) {
		if ($t[6+10*($k-1)] =~ m/^\.*$/) {
			$exon1[6+$k-1] = [];
#			my $gid = '.';
#			if (!exists $gids[$k]{$gid}) { $gids[$k]{$gid} = 0; }
		} else {
# @a: (beg0,end0,name1,beg1,end1,mapped_exons,0,strand1,chain_id,mapped_length)
			my @a = @t[(6+10*($k-1))..(6+10*$k-1)];
			my ($beg,$end) = @a[3..4];
			my @map_exons;
			if ($a[5]!~m/^\.*$/) {
				foreach my $map_pos (split /;/, $a[5]) {
					my @ee = split /,/, $map_pos;
					push @map_exons, \@ee;
					if ($ee[3]=~m/^exon/) { $is_exon = 2; }
					if ($is_exon && $ee[3]=~m/^exon/) {
						my $gid;
						if ($ee[3] =~m/^\.*$/) {
							$gid = '.';
						} else {
							my @ee4 = split /:/, $ee[3];
							$gid = $ee4[5];
						}
						my $b = $ee[1]>$beg?$ee[1]:$beg;
						my $e = $ee[2]<$end?$ee[2]:$end;
						if ($e>$b) { $gids[$k]{$gid} += $e-$b; }
					}
				}
			}
			undef $a[5];
			$a[5] = \@map_exons;
			$exon1[6+$k-1] = \@a;
			$lens[$k] = $a[9];
		}
	}
	$exon1[11] = $is_exon;
	$exon[0]{$ids[0]} = \@exon1;
#	my $map_len = $lens[1]>$lens[2] ? $lens[1] : $lens[2];
#	$gene[0]{$gid0}{map_len} += $map_len;
	if ($gid0 ne '.') { 
		$gene[0]{$gid0}{strand} = $t[5];
		if (!exists $gene[0]{$gid0}{pairs}) { $gene[0]{$gid0}{pairs} = []; }
		if (@{$exon1[6]}>0 && ($exon1[6][1]-$exon1[6][0]>=50 || $exon1[6][1]-$exon1[6][0]>=0.5*($exon1[2]-$exon1[1]))
			&& @{$exon1[7]}>0 && ($exon1[7][1]-$exon1[7][0]>=50 || $exon1[7][1]-$exon1[7][0]>=0.5*($exon1[2]-$exon1[1]))) {
			foreach my $gid1 (sort keys(%{$gids[1]})) {
			foreach my $gid2 (sort keys(%{$gids[2]})) {
				my $l = $gids[1]{$gid1}>$gids[2]{$gid2} ? $gids[1]{$gid1} : $gids[2]{$gid2};
				if ($l<50 && $l<($t[7]-$t[6])*0.2) { next; }
				if ($l>$t[2]-$t[1]) { $l = $t[2]-$t[1]; }
				push @{$gene[0]{$gid0}{pairs}}, [$gid1,$gids[1]{$gid1}, $gid2, $gids[2]{$gid2}, $l]; # l is the matched length on mapped position
			}
			}
		}
	}
}
close IN;

foreach my $eid (keys(%{$exon[0]})) {
	my $e = $exon[0]{$eid};
	if ($e->[3]=~m/exon/) {
		my @ids = split /:/, $e->[3];
		my $gid = $ids[5];
		my $g = $gene[0]{$gid};
		if (!exists $g->{chr}) { $g->{chr} = $e->[0]; }
		if (!exists $g->{strand}) { $g->{strand} = $e->[5]; }
		if (!exists $g->{begin} || $e->[1]<$g->{begin}) {
			$g->{begin} = $e->[1];
		}
		if (!exists $g->{end} || $e->[2]>$g->{end}) {
			$g->{end} = $e->[2];
		}
	}
}


# keep unique pairs
$n=keys(%{$gene[0]});
my $n_good=0;
foreach my $gid0 (sort keys(%{$gene[0]})) {
	my $pairs = $gene[0]{$gid0}{pairs};
	my $i0=0;
	for (my $i=0; $i<@$pairs; $i++) {
		if ($pairs->[$i][0]=~m/^\.*$/ || $pairs->[$i][2]=~m/^\.*$/) { next; }
		if ($i0!=$i) { $pairs->[$i0] = $pairs->[$i];}
		$i0++;
	}

	$#$pairs = $i0-1;
	if (@$pairs==0) { $gene[0]{$gid0}{good}=0; next; }
	if (@$pairs==1) { $gene[0]{$gid0}{good}=1; $n_good++; next; }
	@$pairs = sort {$a->[0] cmp $b->[0] or $a->[2] cmp $b->[2]} @$pairs;

	$i0=0;
	for (my $i=1; $i<@$pairs; $i++) {
		if ($pairs->[$i][0] ne $pairs->[$i0][0] || $pairs->[$i][2] ne $pairs->[$i0][2]) {
			++$i0;
			if ($i0!=$i) { $pairs->[$i0] = $pairs->[$i];}
		} else {
			$pairs->[$i0][1]+=$pairs->[$i][1];
			$pairs->[$i0][3]+=$pairs->[$i][3];
			$pairs->[$i0][4]+=$pairs->[$i][4];
		}

	}
	$#$pairs = $i0;
	@$pairs = sort { -($a->[4]<=>$b->[4]) } @$pairs;
	if (@$pairs==1 || $pairs->[1][4]<0.2*$pairs->[0][4]) { $gene[0]{$gid0}{good}=1; $n_good++; }
	else { $gene[0]{$gid0}{good}=0; next; }
}

#print STDERR "$n_good\n";

# keep unique pairs for sp2
foreach my $gid0 (sort keys(%{$gene[0]})) {
	my $pairs = $gene[0]{$gid0}{pairs};
	if (@$pairs==0) { next; }
	my ($gid1,$l1,$gid2,$l2,$l) = @{$pairs->[0]}[0..4];
	push @{$gene[1]{$gid1}{sp1_match}}, [$gid0, $l1];
	$gene[1]{$gid1}{sp1_match_len} += $l1;
	push @{$gene[1]{$gid2}{sp1_match}}, [$gid0, $l2];
	$gene[1]{$gid2}{sp1_match_len} += $l2;
}
foreach my $gid (sort keys(%{$gene[1]})) {
	my $g = $gene[1]{$gid};
	@{$g->{sp1_match}} = sort {-($a->[1]<=>$b->[1])} @{$g->{sp1_match}};
	my $gid0 = $g->{sp1_match}[0][0];
	my $g0   = $gene[0]{$gid0};
	my $l0   = $g->{sp1_match}[0][1];
	my $i0=0;
	my $s0 = $g0->{strand};
	my $beg0 = $g0->{begin};
	my $end0 = $g0->{end};
	for (my $i=1;$i<@{$g->{sp1_match}}; $i++) {
		my $gid2 = $g->{sp1_match}[$i][0];
		my $g2 = $gene[0]{$gid2};
		if ($g->{sp1_match}[$i][1]<0.2*$l0) {
			for (my $j=$i; $j<@{$g->{sp1_match}}; $j++) {
				my $gid3 = $g->{sp1_match}[$j][0];
				if ($gene[0]{$gid3}{good}!=0) { $gene[0]{$gid3}{good} = 0; $n_good--; }
			}
			$#{$g->{sp1_match}} = $i-1;
			last;
		}
	}
	if (@{$g->{sp1_match}}>1) {
		for (my $i=1;$i<@{$g->{sp1_match}}; $i++) {
			my $gid2 = $g->{sp1_match}[$i][0];
			my $g2 = $gene[0]{$gid2};
			my $beg2 = $g2->{begin};
			my $end2 = $g2->{end};
			my $b = $beg0 > $beg2 ? $beg0 : $beg2;
			my $e = $end0 < $end2 ? $end0 : $end2;
			if ($e-$b<=0) {
				$i0++;
				if ($i!=$i0) { $g->{sp1_match}[$i0] = $g->{sp1_match}[$i]; }
			} else {
				if ($gene[0]{$gid2}{good}!=0) { $gene[0]{$gid2}{good} = 0; $n_good--; }
			}
		}
		$#{$g->{sp1_match}} = $i0;
	}
}
#print STDERR "$n_good\n";

my @good_gid;
my @good_eid;
my %pair;
my @triplets;
foreach my $gid0 (sort keys(%{$gene[0]})) {
	my $pairs = $gene[0]{$gid0}{pairs};
	if ($gene[0]{$gid0}{good}==0) { next; }
	if (@$pairs==0) { next; }
	my ($gid1,$l1,$gid2,$l2,$l) = @{$pairs->[0]}[0..4];
	if (@{$gene[1]{$gid1}{sp1_match}}!=1 || @{$gene[1]{$gid2}{sp1_match}}!=1) {
		$n_good--;
		$gene[0]{$gid0}{good}=0;
		next;
	}
}
#print STDERR "$n_good\n";

foreach my $eid (keys(%{$exon[0]})) {
	my $e = $exon[0]{$eid};
	if ($e->[3]=~m/exon/) {
		my @ids = split /:/, $e->[3];
		my $gid = $ids[5];
		if ($gene[0]{$gid}{good}==1) { push @{$good_eid[0]}, $eid; }
	} else {
		push @{$good_eid[0]}, $eid;
	}
}

foreach my $gid0 (sort keys(%{$gene[0]})) {
	my $pairs = $gene[0]{$gid0}{pairs};
	if ($gene[0]{$gid0}{good}==1) {
		push @{$good_gid[0]}, $gid0;
		my ($gid1,$l1,$gid2,$l2,$l) = @{$pairs->[0]}[0..4];
		$gene[1]{$gid1}{good}=1;
		$gene[1]{$gid2}{good}=1;
		push @triplets, { gids=>[ $gid0, $pairs->[0][0], $pairs->[0][2] ] };
		$gene[0]{$gid0}{triplet_idx} = $#triplets;
		$gene[0]{$gid0}{triplet_idx2} = 0;
		$gene[1]{$gid1}{triplet_idx} = $#triplets;
		$gene[1]{$gid2}{triplet_idx} = $#triplets;
		$gene[1]{$gid1}{triplet_idx2} = 1;
		$gene[1]{$gid2}{triplet_idx2} = 2;
	}
}

my $n0  = @{$good_gid[0]};
my $ne0 = @{$good_eid[0]};

print STDERR "Keep $n0 Gene Triplets, $ne0 exons in species1", "\n";

##########################################
# load sp2 -- sp1 
##########################################
# {{{
open IN, "<$sp2_file" or die;
while(<IN>) {
	if (m/^\s*$/) { next; }
	s/\s+$//;
	my @t = split /\t/;
	my @ids;
	my @begs;
	my @ends;
	$begs[0] = $t[1];
	$ends[0] = $t[2];
	$ids[0] = $t[3];
	$begs[1] = $t[9];
	$ends[1] = $t[10];
	$ids[1] = $t[11];
	my @idss;
	my @gids;
	my $gid0 = '.';
	my @ee = split /:/, $ids[0];
	my $is_exon=0;
	if ($ids[0]=~m/exon/i) { $is_exon = 1; $gid0 = $ee[5]; }

	my $gid1 = '.';
	my $l1=0;
	my @lens = (0,0);
	$lens[0] = $t[2]-$t[1];
	my @exon1 = @t[0..5];
	my $k = 1;
	if ($t[6+10*($k-1)] =~ m/^\.*$/) {
		$exon1[6] = [];
#		my $gid = '.';
#		if (!exists $gids[$k]{$gid}) { $gids[$k]{$gid} = 0; }
	} else {
		my @a = @t[(6+10*($k-1))..(6+10*$k-1)];
		my ($beg,$end) = @a[3..4];
		my @map_exons;
		if ($a[5]!~m/^\.+$/) {
			foreach my $map_pos (split /;/, $a[5]) {
				my @ee = split /,/, $map_pos;
				push @map_exons, \@ee;
				if ($ee[3]=~m/exon/) { $is_exon=2; }
#			my $gid;
#			if ($is_exon && $ee[3]=~m/^exon/) {
#				if ($ee[3] =~m/^\.*$/) {
#					$gid = '.';
#				} else {
#					my @ee4 = split /:/, $ee[3];
#					$gid = $ee4[5];
#				}
#			}
#			my $b = $ee[1]>$beg?$ee[1]:$beg;
#			my $e = $ee[2]<$end?$ee[2]:$end;
#			if ($e>$b) { $gids[$k]{$gid} += $e-$b; }
			}
		}
		undef $a[5];
		$a[5] = \@map_exons;
		$exon1[6] = \@a;
	}
	$exon1[11] = $is_exon;
	$exon[1]{$ids[0]} = \@exon1;
}
close IN;
# }}}

##########################################
# Load sp2 -- sp2 
##########################################
# {{{
open IN, "<$sp22_file" or die;
while(<IN>) {
	if (m/^\s*$/) { next; }
	s/\s+$//;
	my @t = split /\t/;
	my @ids;
	my @begs;
	my @ends;
	$begs[0] = $t[1];
	$ends[0] = $t[2];
	$ids[0] = $t[3];
	$begs[1] = $t[9];
	$ends[1] = $t[10];
	$ids[1] = $t[11];
	my @idss;
	my @gids;
	my $gid0 = '.';
	my @ee = split /:/, $ids[0];
	my $is_exon=0;
	if ($ids[0]=~m/exon/i) { $is_exon = 1; $gid0 = $ee[5]; }

	my $gid1 = '.';
	my $l1=0;
	my @lens = (0,0);
	$lens[0] = $t[2]-$t[1];
	my $k = 1;
	my $exon1 = $exon[1]{$ids[0]};
	if ($t[6+10*($k-1)] =~ m/^\.*$/) {
		$exon1->[7] = [];
#		my $gid = '.';
#		if (!exists $gids[$k]{$gid}) { $gids[$k]{$gid} = 0; }
	} else {
		my @a = @t[(6+10*($k-1))..(6+10*$k-1)];
		my ($beg,$end) = @a[3..4];
		my @map_exons;
		if ($a[5]!~m/^\.+$/) {
			foreach my $map_pos (split /;/, $a[5]) {
				my @ee = split /,/, $map_pos;
				push @map_exons, \@ee;
				if ($ee[3]=~m/^exon/) { $exon1->[11] = 2; }
#			my $gid;
#			if ($is_exon && $ee[3]=~m/^exon/) {
#				if ($ee[3] =~m/^\.*$/) {
#					$gid = '.';
#				} else {
#					my @ee4 = split /:/, $ee[3];
#					$gid = $ee4[5];
#				}
#				my $b = $ee[1]>$beg?$ee[1]:$beg;
#				my $e = $ee[2]<$end?$ee[2]:$end;
#				if ($e>$b) { $gids[$k]{$gid} += $e-$b; }
#			}
			}
		}
		undef $a[5];
		$a[5] = \@map_exons;
		$exon1->[7] = \@a;
#		$lens[$k] = $a[9];
	}
}
close IN;
# }}}

# set gene.exon_len, gene.map_len1,2
foreach my $k (0..1) {
foreach my $eid (keys(%{$exon[$k]})) {
	my $e = $exon[$k]{$eid};
	if ($e->[3]!~m/exon/) { next; }
	my @ids = split /:/, $e->[3];
	my $gid = $ids[5];
	if (!exists $gene[$k]{$gid}) { $gene[$k]{$gid}{good} = 0; }
	$gene[$k]{$gid}{exon_len} += $e->[2]-$e->[1];
	$gene[$k]{$gid}{exon_num} ++;
	if (defined $e->[6] && @{$e->[6]}>0) {
		$gene[$k]{$gid}{map_len1} += $e->[6][1]-$e->[6][0];
	}
	if (defined $e->[7] && @{$e->[7]}>0) {
		$gene[$k]{$gid}{map_len2} += $e->[7][1]-$e->[7][0];
	}
}
}

#######################################################
# assign CNE to nearest genes
#######################################################
# {{{
foreach my $k (0..1) {
foreach my $eid (keys(%{$exon[$k]})) {
	my $e = $exon[$k]{$eid};
	if ($e->[3]=~m/exon/) {
		my @ids = split /:/, $e->[3];
		my $gid = $ids[5];
		my $g = $gene[$k]{$gid};
		if (!exists $g->{chr}) { $g->{chr} = $e->[0]; }
		if (!exists $g->{strand}) { $g->{strand} = $e->[5]; }
		if (!exists $g->{begin} || $e->[1]<$g->{begin}) {
			$g->{begin} = $e->[1];
		}
		if (!exists $g->{end} || $e->[2]>$g->{end}) {
			$g->{end} = $e->[2];
		}
	}
}
}


my @unassign_cne = 0;
my @assigned_cne = 0;
foreach my $k (0..1) {
	my %eids;
	foreach my $eid (keys(%{$exon[$k]})) {
		my $e = $exon[$k]{$eid};
		my $chr = $e->[0];
		push @{$eids{$chr}}, $eid;
	}
	foreach my $chr (keys(%eids)) {
		my $eids1 = $eids{$chr};
		@$eids1 = sort {$exon[$k]{$a}[1]<=>$exon[$k]{$b}[1]} @$eids1;
		my $gid1;
		my $gid2;
		for (my $ei=0; $ei<@$eids1; $ei++) {
			my $eid = $eids1->[$ei];
			my $e = $exon[$k]{$eid};
			if ($e->[3]=~m/exon/) {
				my @ids = split /:/, $e->[3];
				my $gid = $ids[5];
				$e->[9] = $gid;
				$gid1 = $gid;
				undef $gid2;
			} else {
				$e->[9] = $gid1;
				if (!defined $gid2) {
					my $ei2=$ei+1;
					while($ei2<@$eids1) {
						my $eid2 = $eids1->[$ei2];
						my $e2 = $exon[$k]{$eid2};
						if ($e2->[3]=~m/exon/) {
							my @ids2 = split /:/, $e2->[3];
							$gid2 = $ids2[5];
							last;
						}
						$ei2++;
					}
				}
				$e->[10] = $gid2;
				if (!$e->[11]) {
					if (!defined $gid1 && !defined $gid2) { $unassign_cne[$k]++; }
					else { $assigned_cne[$k]++; }
				}
			}
		}
	}
}
# }}}

# count number of lost exons for gene triplets
my @feats = qw(exon CNE CNE5k CNE10k);
my @loss_types = qw(000 001 010 011 100 101 110 111);
foreach my $t (@triplets) {
	foreach my $feat (@feats) {
		foreach my $type (@loss_types) {
			$t->{loss}{$feat}{$type} = [0,0];
		}
	}
}

# count number of lost exons
#foreach my $k (0..1) {
#foreach my $gid (keys(%{$gene[$k]})) {
#	foreach my $feat (qw(exon CNE CNE5k)) {
#		$gene[$k]{$gid}{"${feat}_lost12"}=0;
#		$gene[$k]{$gid}{"${feat}_lost1"}=0;
#		$gene[$k]{$gid}{"${feat}_lost2"}=0;
#		$gene[$k]{$gid}{"${feat}_lost0"}=0;
#	}
#}
#}

my @intron_cne_count; # in intron
my @distant_cne_count;
my @cne5k_count;
my @cne10k_count;
my @pass_eid;
foreach my $k (0..1) {
foreach my $eid (keys(%{$exon[$k]})) {
	if (exists $pass_eid[$k]{$eid}) { next; }
	$pass_eid[$k]{$eid}++;
	my $e = $exon[$k]{$eid};
	my $loss1=0;
	my $loss2=0;
	if (!defined $e->[6] || @{$e->[6]}==0) { $loss1=1; }
	else {
		my $map_span = $e->[6][1]-$e->[6][0];
		if ($map_span < 50 && $map_span < 0.5*($e->[2]-$e->[1])) { $loss1=1; }
	}
	if (!defined $e->[7] || @{$e->[7]}==0) { $loss2=1; }
	else {
		my $map_span = $e->[7][1]-$e->[7][0];
		if ($map_span < 50 && $map_span < 0.5*($e->[2]-$e->[1]))
		{ $loss2=1; }
	}

	if (@{$e->[6]}>0) {
		foreach my $ee (@{$e->[6][5]}) {
			my $map_eid = $ee->[3];
			$pass_eid[1-$k]{$map_eid}++;
		}
	}
	if (@{$e->[7]}>0) {
		foreach my $ee (@{$e->[7][5]}) {
			my $map_eid = $ee->[3];
			$pass_eid[1-$k]{$map_eid}++;
		}
	}

	my $gid;
	my @feats;
	if ($e->[3]=~m/exon/) {
		push @feats, 'exon';
		my @ids = split /:/, $e->[3];
		$gid = $ids[5];
	} elsif ($e->[11]==2) {
		push @feats, 'exon';
		if (defined $e->[9]) {
			if (defined $e->[10]) {
				if ($e->[9] eq $e->[10]) { $gid = $e->[9]; }
				else {
					my $gid1 = $e->[9];
					my $gid2 = $e->[10];
					my ($d1,$d2);
					$d1 = $e->[1]-$gene[$k]{$gid1}{end};
					$d2 = $gene[$k]{$gid2}{begin}-$e->[2];
					if ($d1*2<$d2 && $d1<=5000) { $gid = $e->[9]; undef $e->[10];}
					elsif ($d2*2<$d1 && $d2<=5000) { $gid=$e->[10]; undef $e->[9]; }
					else { next; }
				}
			} else {
				$gid = $e->[9];
			}
		} else {
			if (defined $e->[10]) {
				$gid = $e->[10];
			} else {
				next;
			}
		}
	} elsif ($e->[11]==0) {
		push @feats, 'CNE';
		my $j = 0;
		my $d = 0;
		if (defined $e->[9]) {
			if (defined $e->[10]) {
				my $gid1 = $e->[9];
				my $gid2 = $e->[10];
				my ($d1,$d2);
				$d1 = $e->[1]-$gene[$k]{$gid1}{end};
				$d2 = $gene[$k]{$gid2}{begin}-$e->[2];
				$j = $d1<$d2 ? 9 : 10;
				if ($gid1 eq $gid2) { $intron_cne_count[$k]++; }
			} else {
				my $gid1 = $e->[9];
				$d = $e->[1]-$gene[$k]{$gid1}{end};
				$j = 9;
			}
		} else {
			if (defined $e->[10]) {
				$j = 10;
				my $gid2 = $e->[10];
				$d = $gene[$k]{$gid2}{begin}-$e->[2];
			} else { next; }
		}
		$gid = $e->[$j];
		if ($d<=5000) { push @feats, 'CNE5k'; $cne5k_count[$k]++; }
		if ($d<=10000) { push @feats, 'CNE10k'; $cne10k_count[$k]++; }
		else { $distant_cne_count[$k]++; }
	}

	my $g = $gene[$k]{$gid};
	if (!exists $g->{good} || $g->{good}==0) { next; }
	my $tidx = $g->{triplet_idx};
	my $tidx2 = $g->{triplet_idx2};
	my $t = $triplets[$tidx];

	my $type;
	if ($k==0) {
		if ($loss1==0) {
			if ($loss2==0) { $type='000'; }
			else { $type='001'; }
		} else {
			if ($loss2==0) { $type='010'; }
			else { $type='011'; }
		}
	} else {
		if ($loss1==0) {
			if ($loss2==0) { $type='000'; }
			else {
				if ($tidx2==1) { $type='001'; }
				elsif ($tidx2==2) { $type='010'; }
			}
		} else {
			if ($loss2==0) {
				$type = '100';
			} else {
				if ($tidx2==1) { $type='101'; }
				elsif ($tidx2==2) { $type='110'; }
			}
		}
	}

	foreach my $feat (@feats) {
		$t->{loss}{$feat}{$type}[0]++;
		$t->{loss}{$feat}{$type}[1]+=$e->[2]-$e->[1];
	}
}
}

#my @lost_num;
#my @lost_cne_num;
#my @lost_cne_5k_num;
#foreach my $k (0..1) {
#foreach my $gid (keys(%{$gene[$k]})) {
#	my $g = $gene[$k]{$gid};
#	if (!$g->{good}) { next; }
#	foreach my $feat (qw(exon CNE CNE5k)) {
#		my $loss;
#		if (exists $g->{"${feat}_lost1"}) {
#			if (exists $g->{"${feat}_lost2"}) {
#				$loss=3;
#			} else {
#				$loss=1;
#			}
#		} else {
#			if (exists $g->{"${feat}_lost2"}) {
#				$loss=2;
#			} else {
#				$loss=0;
#			}
#		}
#		$lost_num[$k]{$feat}[$loss]++;
#	}
#}
#}
#
#foreach my $k (0..1) {
#foreach my $i (0..3) {
#	foreach my $feat (qw(exon CNE CNE5k)) {
#		if (! defined $lost_num[$k]{$feat}[$i]) {$lost_num[$k]{$feat}[$i]=0;}
#		print "${feat}_loss\t$k\t$i\t$lost_num[$k]{$feat}[$i]\n";
#	}
#}
#}

my $file = "${out_prefix}triplet.stat.txt";
if (-f $file) { die "File $file exists.\n"; }
open OUT, ">$file";
print OUT "#ID1\tID2a\tID2b";
foreach my $feat (@feats) {
	foreach my $type (@loss_types) {
		print OUT "\t${feat}_${type}_num";
	}
	foreach my $type (@loss_types) {
		print OUT "\t${feat}_${type}_len";
	}
}
print OUT "\n";
$k = 0;
foreach my $t (@triplets) {
	print OUT join("\t", @{$t->{gids}});
	foreach my $feat (@feats) {
		foreach my $type (@loss_types) {
			print OUT "\t", $t->{loss}{$feat}{$type}[0];
		}
		foreach my $type (@loss_types) {
			print OUT "\t", $t->{loss}{$feat}{$type}[1];
		}
	}
	print OUT "\n";
}
close OUT;
#if ($copy==1) {
#	foreach my $gid0 (sort keys(%gene)) {
#	foreach my $gid1 (sort keys(%{$gene{$gid0}})) {
#		print join("\t", ($gid0,$gid1, $gene{$gid0}{$gid1})), "\n";
#	}
#	}
#} else {
#	foreach my $gid0 (sort keys(%gene)) {
#	foreach my $gid1 (sort keys(%{$gene{$gid0}})) {
#	foreach my $gid2 (sort keys(%{$gene{$gid0}{$gid1}})) {
#		my $a = $gene{$gid0}{$gid1}{$gid2};
#		print join("\t", ($gid0,$gid1, $a->[0], $gid2, $a->[1])), "\n";
#	}
#	}
#	}
#}

