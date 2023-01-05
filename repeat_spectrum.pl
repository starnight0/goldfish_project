#!/usr/bin/perl -I$HOME/my_program3/src/perl_pm
use strict;
use warnings;
use Data::Dumper;

my $col = 4;
$col--;
my $sep = '__'; # species, family, class
my %repeat_count;
my %repeat_bp;
my %repeat_count2;
my %repeat_bp2;
while(<STDIN>) {
	if (m/^#/) { next; }
	if (m/^\s*$/) { next; }
	my @tab = split /\t/;
	if ($tab[$col]=~m/\?/)  { next; } # ignore uncertained classification
	$tab[$col]=~s/^\s+//;
	$tab[$col]=~s/\s+$//;
	my @t = split $sep, $tab[$col];
	if (@t==2) {
		$repeat_count{$t[1]}{ALL}{$t[0]}++;
		$repeat_bp{$t[1]}{ALL}{$t[0]}+=$tab[2]-$tab[1];
		$repeat_count2{$t[1]}{ALL}++;
		$repeat_bp2{$t[1]}{ALL}+=$tab[2]-$tab[1];
	} elsif (@t>=3) {
		$repeat_count{$t[1]}{$t[2]}{$t[0]}++;
		$repeat_bp{$t[1]}{$t[2]}{$t[0]} += $tab[2]-$tab[1];
		$repeat_count2{$t[1]}{$t[2]}++;
		$repeat_bp2{$t[1]}{$t[2]}+=$tab[2]-$tab[1];
	}
}
print join("\t", qw(D2 D1 D0 Count_210 bp_210 Count_21 bp_21)), "\n";
foreach my $s0 (sort keys(%repeat_count)) {
	foreach my $s1 (sort keys(%{$repeat_count{$s0}})) {
		foreach my $s2 (sort keys(%{$repeat_count{$s0}{$s1}})) {
			my $n=$repeat_count{$s0}{$s1}{$s2};
			my $m=$repeat_bp{$s0}{$s1}{$s2};
			my $n2=$repeat_count2{$s0}{$s1};
			my $m2=$repeat_bp2{$s0}{$s1};
			print join("\t", ($s0,$s1,$s2, $n,$m, $n2,$m2)), "\n";
		}
	}
}
