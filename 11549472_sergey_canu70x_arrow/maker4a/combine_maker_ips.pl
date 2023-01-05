#!/usr/bin/perl -I$HOME/my_program3/src/perl_pm
use strict;
use warnings;
use Data::Dumper;
use czl_io::base_io;

sub usage()
{
print<<EOF;
Usage:
  combine_maker_ips.pl MAKER_GFF IPS_TSV OUT
EOF
}

if ($#ARGV<2) {
	usage();
	exit 0;
}

use constant {
	IPS_ID        => 0,
	IPS_MD5       => 1,
	IPS_SEQLEN    => 2,
	IPS_ANALYSIS  => 3,
	IPS_ACC       => 4,
	IPS_DESC      => 5,
	IPS_START     => 6,
	IPS_STOP      => 7,
	IPS_SCORE     => 8,
	IPS_STATUS    => 9,
	IPS_DATE      => 10,
	IPS_INTERPRO_ACC   => 11,
	IPS_INTERPRO_DESC  => 12,
	IPS_GO        => 13,
	IPS_PATHWAY   => 14
};

my $igff    = $ARGV[0];
my $ips_tsv = $ARGV[1];
my $ogff    = $ARGV[2];

my $fin;
my $fout;

my %ips;
$fin = base_io::czl_open($ips_tsv, "r");
while (<$fin>) {
	if (m/^#/) { next; }
	s/^\s*//; s/\s*$//;
	my @tab = split "\t";
	my $rna = $tab[0];
	if ($#tab>=IPS_INTERPRO_ACC && $tab[IPS_INTERPRO_ACC]=~m/\S/) {
		my $acc = $tab[IPS_INTERPRO_ACC];
		if (exists $ips{$rna} && exists $ips{$rna}{interpro}{$acc}) {
		} else {
			$ips{$rna}{interpro}{$acc}{print} = "$acc:$tab[IPS_INTERPRO_DESC]";
			$ips{$rna}{interpro}{$acc}{print} =~s/\s/_/g;
		}
		$ips{$rna}{interpro}{$acc}{count}++;
	}

	if ($#tab>=IPS_GO && $tab[IPS_GO]=~m/\S/) {
		foreach my $go (split /\|/, $tab[IPS_GO]) {
			if (exists $ips{$rna} && exists $ips{$rna}{go}{$go}) {
			} else {
				$ips{$rna}{go}{$go}{print} = $go;
			}
			$ips{$rna}{go}{$go}{count}++;
		}
	}

	if ($#tab>=IPS_PATHWAY && $tab[IPS_PATHWAY]=~m/\S/) {
		foreach my $aa (split /\|/, $tab[IPS_PATHWAY]) {
			my ($ana,$acc) = split ":", $aa;
			if (exists $ips{$rna} && exists $ips{$rna}{$ana}{$acc}) {
			} else {
				$ips{$rna}{$ana}{$acc}{print} = $acc;
			}
			$ips{$rna}{$ana}{$acc}{count}++;
		}
	}

	my $ana = $tab[IPS_ANALYSIS];
	my $acc = $tab[IPS_ACC];
	if (!exists $ips{$rna}{$ana}{$acc}) {
		my $desc = $tab[IPS_DESC];
		$desc =~ s/^\s*//;
		$desc =~ s/\s*$//;
		$desc =~ s/\s/_/g;
		if (exists $ips{$rna}{$ana} && exists $ips{$rna}{$ana}{$acc}) {
		} else {
			if ($desc=~/^\s*/) {
				$ips{$rna}{$ana}{$acc}{print} = "$acc";
			} else {
				$ips{$rna}{$ana}{$acc}{print} = "$acc:$desc";
			}
		}
		$ips{$rna}{$ana}{$acc}{count}++;
	}
}
close $fin;

$fin =base_io::czl_open($igff, "r");
$fout=base_io::czl_open($ogff, "w");
while(<$fin>) {
	if (m/^#/) { print $fout; next; }
	chomp;
	my @tab = split "\t";
	my $id0;
	foreach my $attr (split /\s*;\s*/, $tab[8]) {
		my ($u,$v) = split /\s*=\s*/, $attr;
		if ($u eq "ID") {
			$id0 = $v;
		}
	}
	if ($tab[2] eq "mRNA") {
		if (exists $ips{$id0}) {
			foreach my $ana (sort keys(%{$ips{$id0}})) {
				my @accs = sort keys(%{$ips{$id0}{$ana}});
				my @prints;
				foreach my $acc (@accs) {
					push @prints, $ips{$id0}{$ana}{$acc}{print};
				}
				$tab[8] .= ";$ana=" . join("|",@prints);
			}
		}
	}
	print $fout join("\t",@tab), "\n";
}
close $fin;
close $fout;
