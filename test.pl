use threads;
use threads::shared;
use strict;
use warnings;

my %a :shared;
%a = ( A=>1,
		B=>2,
		);
my $b :shared;

my $th = threads->create(sub { $a{C}=3; $a{D} = shared_clone([4,5]); print $a{D}, "\n"; $b=\%a; print $b->{D}, "\n"; return \%a;});
$th->join();

$b->{D}[0]=6;
print \%a, "\n";
foreach my $keys (sort keys(%a)) {
	print "$keys\t$a{$keys}\n";
}
print join(',',@{$a{D}}), "\n";
print "\n";
print "\n";
print $b, "\n";
foreach my $keys (sort keys(%$b)) {
	print "$keys\t$b->{$keys}\n";
}
print join(',',@{$b->{D}}), "\n";

print "\n";
print "\n";
my %c=(A=>1,B=>2);
$c{E}=[4,5];
my $d=\%c;
print \%c, "\n";
foreach my $keys (sort keys(%c)) {
	print "$keys\t$c{$keys}\n";
}

print $d, "\n";
foreach my $keys (sort keys(%$d)) {
	print "$keys\t$d->{$keys}\n";
}

