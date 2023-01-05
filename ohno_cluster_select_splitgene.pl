#!/usr/bin/perl -I$HOME/my_program3/src/perl_pm

use strict;
use warnings;
use czl_io::base_io;

my $anchor_file = $ARGV[0];
my $clust_file=$ARGV[1];
# my $align_file = $ARGV[2];

my $fin;

my %anchor;
my @anchors;
$fin=base_io::czl_open($anchor_file, "r");
while(<$fin>) {
	chomp();
	my @t=split /\t/,$_,-1;
	if (@t>0) {
		$anchor{$t[0]}{$t[8]} = \@t;
		push @anchors, \@t;
	}
}
close $fin;

#my %align;
#$fin=base_io::czl_open($align_file, "r");
#while(<$fin>) {
#	chomp();
#	my @t=split /\t/,$_,-1;
#	$align{$t[0]}{$t[1]} = \@t;
#}
#close $fin;


$fin=base_io::czl_open($clust_file, "r");
while(<$fin>) {
	chomp;
	my @t = split /\t/;
	foreach my $a (split /;/, $t[3]) {
		foreach my $a2 (split /,/, $a) {
			if ($a2=~m/^\((.*)\)$/) {
				my @g = split /\|/, $1;
				my ($sp0, $chr0, $gid0, $gname0) = split /:/, $g[0];
				for (my $i=1; $i<@g; $i++) {
					my $g1 = $g[$i];
					my ($sp1, $chr1, $gid1, $gname1) = split /:/, $g1;
					$anchor{$sp1}{$gid1}[10] = $gid0; # fix 'to_dup' in anchor file
				}
			}
		}
	}
}
close $fin;

foreach my $a (@anchors) {
	print join("\t", @$a), "\n";
}
