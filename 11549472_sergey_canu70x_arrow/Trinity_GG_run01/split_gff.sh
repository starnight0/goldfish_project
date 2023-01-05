cat  all.f.ex2.ncrna.cuffcompare.combined.gtf | perl -e 'while(<>) {
if (m/^#/ || m/^\s*$/) { next; } chomp;
my @t=split "\t"; my $chr = $t[0];

}'
