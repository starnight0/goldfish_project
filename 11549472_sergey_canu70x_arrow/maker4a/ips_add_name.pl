#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use czl_io::base_io;
use czl_annot qw(load_kegg);

sub usage
{
print <<EOF;
Usage:
  $0 IPS_TSV OUT
EOF
}

if ($#ARGV<1) {
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
# interprotscan output format:
# Protein Accession (e.g. P51587)
# Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
# Sequence Length (e.g. 3418)
# Analysis (e.g. Pfam / PRINTS / Gene3D)
# Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
# Signature Description (e.g. BRCA2 repeat profile)
# Start location
# Stop location
# Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
# Status - is the status of the match (T: true)
# Date - is the date of the run
# (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
# (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
# (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
# (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
#
# Analysis can be:
# CDD               have desc
# Coils             no desc
# Gene3D            no desc
# Hamap             have desc
# MobiDBLite        have desc
# PANTHER           no desc
# PIRSF             no desc
# PRINTS            have desc
# Pfam              have desc
# ProSitePatterns   have desc
# ProSiteProfiles   have desc
# SFLD              have desc
# SMART             no desc
# SUPERFAMILY       no desc
# TIGRFAM           no desc
my $ips_tsv = $ARGV[0];
my $out_tsv = $ARGV[1];
my $fin;
my $fout;

my $kegg_dir = "/home/chenz11/data/KEGG";
my $kegg = czl_annot::load_kegg($kegg_dir);
# File:  ReactomePathwarys.txt
# ID,desc,species, separated by TAB
my $reactome_file="/home/chenz11/data/reactome/ReactomePathways.txt";

# File:  SequenceAssociationPathway3.5.txt
# Pathway accession
# Pathway name
# Pathway component accession
# Pathway component name
# UniProt ID
# Protein definition
# Confidence code
# Evidence
# Evidence type (e.g., PubMed, OMIM)
# PANTHER subfamily ID (the SF to which the sequence belongs to)
# PANTHER subfamily name
my $panther_file="/home/chenz11/data/pantherdb/pathway/3.5/SequenceAssociationPathway3.5.txt";

# File: go.obo
my $go_obo="/home/chenz11/data/GO/go.obo";

my $go = czl_annot::load_go_obo($go_obo);

my %reactome;
$fin = base_io::czl_open($reactome_file, "r") or die "Fail to opern $reactome_file\n";
while (<$fin>) {
	my ($id, $name, $sp) = split "\t";
	$reactome{$id} = { name=>$name, species=>$sp};
}
close $fin;

my %panther;
$fin = base_io::czl_open($panther_file, "r");
while (<$fin>) {
	chomp;
	my ($pw_id, $pw_name, $pw_subid, $pw_subname, $uniprot_id, $uniprot_name, $evd, $evd_type, $sf_id, $sf_name) = split "\t";
	if ( !exists $panther{pathway}{$pw_id} ) {
		$panther{pathway}{$pw_id} = { name=>$pw_name};
	}
	if ( !exists $panther{pathway}{$pw_subid} ) {
		$panther{pathway}{$pw_subid} = { name=>$pw_subname, parent=>$pw_id};
	}
	if ( !exists $panther{pathway}{$pw_subid} ) {
		$panther{family}{$sf_id} = { name=>$sf_name }
	}
	if (!exists $panther{family}{$sf_id}{pathway}{$pw_id}) {
		$panther{family}{$sf_id}{pathway}{$pw_id}++;
		push @{$panther{pathway}{$pw_id}{family}}, $sf_id;
	}
	if (!exists $panther{family}{$sf_id}{pathway}{$pw_subid}) {
		$panther{family}{$sf_id}{pathway}{$pw_subid}++;
		push @{$panther{pathway}{$pw_subid}{family}}, $sf_id;
	}
}
close $fin;


$fin = base_io::czl_open($ips_tsv, "r") or die "Fail to opern $ips_tsv\n";
$fout= base_io::czl_open($out_tsv, "w");
while (<$fin>) {
	if (m/^#/) { print $fout $_; next; }
	chomp;
	my @tab = split "\t";
	my $ana = $tab[IPS_ANALYSIS];
	if (defined $tab[IPS_GO] && $tab[IPS_GO]=~/\S/) {
		my @goterms = split /\|/, $tab[IPS_GO];
		# START add name and domain for goterms
		foreach my $goterm (@goterms) {
			my $d;
			my $domain;
			foreach my $d1 (keys(%$go)) {
				if (exists $go->{$d1}{$goterm}) { $domain=$d1; last; }
			}
			my $go_name = "";
			if (defined $domain) {
				$d = $domain;
				$d =~ s/^(.)[A-Za-z]+_(.).*$/$1$2/;
				if (exists $go->{$domain}{$goterm} && exists $go->{$domain}{$goterm}{name} ) {
					$go_name = $go->{$domain}{$goterm}{name};
					$go_name =~ s/\s/_/g;
				}
			} else {
				$d="";
			}
			$goterm .= ":$d:$go_name";
		}
		$tab[IPS_GO] = join '|', @goterms;
	}
	# END
	
	if (defined $tab[IPS_PATHWAY] && $tab[IPS_PATHWAY]=~/\S/) {
		my @pathways = split /\|/, $tab[IPS_PATHWAY];
		foreach my $p (@pathways) {
			my ($u,$v) = split /\s*:\s*/, $p;
			if ($u =~ m/^Reactome/i) {
				my $name = "";
				if (exists $reactome{$v}) {
					$name=$reactome{$v}{name};
					$name =~ s/\s/_/g;
				}
				$p = "$u:$v:$name";
			} elsif ($u =~ m/^KEGG/i) {
				my ($id,$name) = split /\+/, $v;
				my $desc="";
				if (exists $kegg->{pathway}{"map:$id"}) {
					$desc = $kegg->{pathway}{"map:$id"}{desc};
					$desc =~ s/\s/_/g;
				}
				$p="$u:$v:$desc";
			} elsif ($u =~ m/^MetaCyc/i) {
			} else {
			}
		}
		$tab[IPS_PATHWAY] = join '|', @pathways;
	}
	print $fout join("\t",@tab), "\n";
}
close $fin;
close $fout;
