use strict;
use warnings;

if ($#ARGV<1) { die; }
open IN1, "<$ARGV[0]" or die;
open IN2, "<$ARGV[1]" or die;

my %to_clust;
while (<IN1>) {
	my @t=split /\t/;
	my $cid = $t[0];
	$t[3]=~s/\(//g;
	$t[3]=~s/\)//g;
	my @genes = split /[,;]/, $t[3];
	foreach my $g (@genes) {
		my ($sp, $chr, $gid, $gname) = split /:/, $g;
		$to_clust{$gid} = $cid;
	}
}
close IN1;


my $n_same=0;
my $n_diff=0;
my $n_nan=0;
while (<IN2>) {
	my @t=split /\t/;
	my $gid1 = $t[0];
	my $diff = 0;
	my $nan=0;
	my $cid;
	if (exists $to_clust{$gid1}) {
		$cid = $to_clust{$gid1};
	} else {
		foreach my $a (@t[1..$#t]) {
			my @t1 = split /,/, $a;
			foreach my $a1 (@t1) {
				my @t2 = split /:/, $a1;
				my $gid2 = $t2[5];
				if (!exists $to_clust{$gid2}) { $nan++; }
				else {
					$cid = $to_clust{$gid2};
					last;
				}
			}
			if (defined $cid) { last; }
		}
	}
	foreach my $a (@t[1..$#t]) {
		my @t1 = split /,/, $a;
		foreach my $a1 (@t1) {
			my @t2 = split /:/, $a1;
			my $gid2 = $t2[5];
			if (exists $to_clust{$gid2}) {
				if ($to_clust{$gid2} != $cid) { $diff++; }
			}
		}
	}
	if ($diff>0) { $n_diff++; }
	else { $n_same++; }
	if ($nan>0) { $n_nan++; }
}
close IN2;

print "same=$n_same\n";
print "diff=$n_diff\n";
print "nan=$n_nan\n";
