use strict;
use warnings;
use Data::Dumper;

my $in_file = $ARGV[0];
my $m6_file = $ARGV[1];
my $out_prefix = $ARGV[2];
my $no_btop=0;
if (defined $ARGV[3] && $ARGV[3] eq 'noBTOP') {$no_btop=1;}
my $file;

if ($in_file=~m/\.gz$/) {
	open IN, "gzip -dc $in_file |";
} else {
	open IN, "<$in_file";
}
my @gids;
my %pairs;
while(<IN>) {
	if (m/^\s*$/) { next; }
	if (m/^#/) { next; }
	s/\s+$//;
	my @t = split /\t/;
	push @gids, \@t;
	my %a;
	for (my $i1=0; $i1<=$#t; $i1++) {
	for (my $i2=0; $i2<=$#t; $i2++) {
		$pairs{$t[$i1]}{$t[$i2]}= [];
	}
	}
}
close IN;


if ($m6_file=~m/\.gz$/) {
	open IN, "gzip -dc $m6_file |";
} else {
	open IN, "<$m6_file";
}
while(<IN>) {
	if (m/^\s*$/) { next; }
	if (m/^#/) { next; }
	s/\s+$//;
	my @t = split /\t/;
	my $id1 = $t[0];
	my $id2 = $t[1];
	$id1=~s/\.[0-9]+$//;
	$id2=~s/\.[0-9]+$//;
	if (exists $pairs{$id1}{$id2}) {
		$pairs{$id1}{$id2} = \@t;
	}
}
close IN;

$file = "${out_prefix}m6";
if (-f $file) { die "File $file exists.\n"; }
open OUT, ">$file";
foreach my $gids1 (@gids) {
	for (my $i1=0; $i1<=$#$gids1; $i1++) {
	for (my $i2=0; $i2<=$#$gids1; $i2++) {
		my $id1 = $gids1->[$i1];
		my $id2 = $gids1->[$i2];
		my $p = $pairs{$id1}{$id2};
		if (@$p>0) {
			print OUT join("\t", @$p), "\n";
		}
	}
	}
}
close OUT;

my $n_id = @{$gids[0]};

$file = "${out_prefix}align_stat.txt";
if (-f $file) { die "File $file exists.\n"; }
open OUT, ">$file";
foreach my $gids1 (@gids) {
	print OUT join("\t", @$gids1);
	for (my $i1=0; $i1<=$#$gids1; $i1++) {
	for (my $i2=0; $i2<=$#$gids1; $i2++) {
		my $id1 = $gids1->[$i1];
		my $id2 = $gids1->[$i2];
		my $p = $pairs{$id1}{$id2};
		my ($size1,$size2,$cov1,$cov2,$iden,$align_len,$bit,$score);
		if (@$p>0) {
			$iden = $p->[2];
			$align_len = $p->[3];
			$bit  = $p->[11];
			if ($no_btop) {
				$size1= $p->[12];
				$size2 = $p->[13];
				$score = $p->[14];
				$cov1 = $p->[17];
				$cov2 = $p->[18];
			} else {
				$size1 = $p->[13];
				$size2 = $p->[14];
				$score = $p->[15];
				$cov1 = $p->[18];
				$cov2 = $p->[19];
			}
		} else {
			$p = $pairs{$id2}{$id1};
			if (@$p>0) {
				$iden = $p->[2];
				$align_len = $p->[3];
				$bit  = $p->[11];
				if ($no_btop) {
					$size1 = $p->[13];
					$size2 = $p->[12];
					$score = $p->[14];
					$cov1 = $p->[18];
					$cov2 = $p->[17];
				} else {
					$size1 = $p->[14];
					$size2 = $p->[13];
					$score = $p->[15];
					$cov1 = $p->[19];
					$cov2 = $p->[18];
				}
			}
		}
		if (@$p>0) {
			print OUT "\t", join("\t", ($size1,$size2,$cov1,$cov2,$iden,$align_len,$bit,$score));
		} else {
			print OUT "\t", join("\t", (0) x 8);
		}
	}
	}
	print OUT "\n";
}
close OUT;
