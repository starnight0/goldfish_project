#!/bin/perl -I$HOME/my_program3/src/perl_pm
use strict;
use warnings;
use Data::Dumper;
use czl_io::gff;
use czl_io::base_io;

sub usage()
{
print<<EOF;
Usage:
  assign_name.pl INPUT_GFF INPUT_GFF_FOR_NAME Tid_Tname_Ttype_Gid_Gname_Gtype.tab OUT_GFF
EOF
}

if (@ARGV<1) { &usage(); exit(0); }

my $igff =$ARGV[0];
my $map_file=$ARGV[1];
my $gt_file=$ARGV[2];
my $ogff =$ARGV[3];

my %annot;
my %gene1;
my $in=base_io::czl_open($gt_file, "r") or die "Fail to open $gt_file\n";
while(<$in>) {
	if (m/^#/) { next; }
	chomp;
	my @tab = split "\t";
	my ($tid,$tname,$ttype,$gid,$gname,$gtype) = @tab[0..5];
	$tid=~s/\.[0-9]+$//;
	$gid=~s/\.[0-9]+$//;
	$annot{$tid} = { gid=>$gid, gname=>$gname, gtype=>$gtype};
	$gene1{$gid} = { name=>$gname, type=>$gtype};
}
close $in;

my %map;
my %gmap;
$in=base_io::czl_open($map_file, "r") or die "Fail to open $map_file\n";
while(<$in>) {
	if (m/^#/) { next; }
	chomp;
	my @tab = split "\t";
	my $rna_id = $tab[0];
	my $gene_id=$rna_id;
	$gene_id=~s/-mRNA.*$//;
	my $id1 = $tab[1];
	$id1=~s/\.[0-9]+$//;
	if (! exists $annot{$id1}) { next; }
	my $gid1 = $annot{$id1}{gid};
	$gid1=~s/\.[0-9]+$//;
	if (exists $map{$rna_id} && exists $map{$rna_id}{$id1}) {
		if ($tab[3]>$map{$rna_id}{$id1}{score}) {
			$map{$rna_id}{$id1}{iden }=$tab[2]; 
			$map{$rna_id}{$id1}{score}=$tab[3];
		}
	} else {
		$map{$rna_id}{$id1}{iden}  = $tab[2];
		$map{$rna_id}{$id1}{score} = $tab[3];
	}
	if (exists $gmap{$gene_id} && exists $gmap{$gene_id}{$gid1}) {
		if ($tab[3]>$gmap{$gene_id}{$gid1}{score}) {
			$gmap{$gene_id}{$gid1}{iden }=$tab[2]; 
			$gmap{$gene_id}{$gid1}{score}=$tab[3];
		}
	} else {
		$gmap{$gene_id}{$gid1}{iden}  = $tab[2];
		$gmap{$gene_id}{$gid1}{score} = $tab[3];
	}
}
close $in;

foreach my $id0 (keys(%map)) {
	foreach my $id1 (keys(%{$map{$id0}})) {
		$map{$id0}{$id1}{name} = "$annot{$id1}{gid}__$annot{$id1}{gname}__$annot{$id1}{gtype}";
	}
}
foreach my $id0 (keys(%map)) {
	my $max_score = 0;
	foreach my $id1 (keys(%{$map{$id0}})) {
		if ($map{$id0}{$id1}{score} > $max_score) { $max_score = $map{$id0}{$id1}{score};}
	}
	foreach my $id1 (keys(%{$map{$id0}})) {
		if ($map{$id0}{$id1}{score} < 0.8*$max_score) { delete $map{$id0}{$id1}; }
	}
}

foreach my $id0 (keys(%gmap)) {
	foreach my $id1 (keys(%{$gmap{$id0}})) {
		$gmap{$id0}{$id1}{name} = "${id1}__$gene1{$id1}{name}__$gene1{$id1}{type}";
	}
}
foreach my $id0 (keys(%gmap)) {
	my $max_score = 0;
	foreach my $id1 (keys(%{$gmap{$id0}})) {
		if ($gmap{$id0}{$id1}{score} > $max_score) { $max_score = $gmap{$id0}{$id1}{score};}
	}
	foreach my $id1 (keys(%{$gmap{$id0}})) {
		if ($gmap{$id0}{$id1}{score} < 0.8*$max_score) { delete $gmap{$id0}{$id1}; }
	}
}

$in=base_io::czl_open($igff, "r");
my $out=base_io::czl_open($ogff, "w");
while(<$in>) {
	if (m/^#/) { next; }
	chomp;
	my @tab = split "\t";
	my $id0;
	foreach my $attr (split /\s*;\s*/, $tab[8]) {
		my ($u,$v) = split /\s*=\s*/, $attr;
		if ($u eq "ID") {
			$id0 = $v;
		}
	}
	if ($tab[2] eq "gene") {
		if (exists $gmap{$id0}) {
			$tab[8] =~ s/([;\t])Name=/$1MakerName=/;
			my $name;
			my $names;
			if (keys(%{$gmap{$id0}})==1) {
				my @id1 = keys(%{$gmap{$id0}});
				$name = $gmap{$id0}{$id1[0]}{name};
				$tab[8] .= ";Name=$name";
			} else {
				foreach my $id1 ( sort {$gmap{$id0}{$b}{score} <=> $gmap{$id0}{$a}{score}} keys(%{$gmap{$id0}})) {
					if (defined $name) { $names .= ","; }
					else { $name = $gmap{$id0}{$id1}{name}; }
					$names .= $gmap{$id0}{$id1}{name};
				}
				$tab[8] .= ";Name=$name,...;names=$names";
			}
		}
		print $out join("\t",@tab), "\n";
	} elsif ($tab[2] eq "mRNA") {
		if (exists $map{$id0}) {
			$tab[8] =~ s/([;\t])Name=/$1MakerName=/;
			my $name;
			my $names;
			if (keys(%{$map{$id0}})==1) {
				my @id1 = keys(%{$map{$id0}});
				$name = "$id1[0]_$map{$id0}{$id1[0]}{name}";
				$tab[8] .= ";Name=$name";
			} else {
				foreach my $id1 ( sort {$map{$id0}{$b}{score} <=> $map{$id0}{$a}{score}} keys(%{$map{$id0}})) {
					if (defined $name) { $names .= ","; }
					else { $name = "${id1}_$map{$id0}{$id1}{name}"; }
					$names .= "${id1}_$map{$id0}{$id1}{name}";
				}
				$tab[8] .= ";Name=$name,...;names=$names";
			}
		}
		print $out join("\t",@tab), "\n";
	} else {
		print $out "$_\n";
	}
}
close $in;
close $out;
