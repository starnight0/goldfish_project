#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
use czl_io::base_io;
use czl_io::fasta;

use constant {
	M6_QID   =>0,
	M6_TID   =>1,
	M6_PIDEN =>2,
	M6_ALEN  =>3,
	M6_MIS   =>4,
	M6_GAPO  =>5,
	M6_QBEGIN=>6,
	M6_QEND  =>7,
	M6_TBEGIN=>8,
	M6_TEND  =>9,
	M6_EVALUE=>10,
	M6_BIT   =>11,
	M6_BTOP  =>12,
	M6_QSIZE =>13,
	M6_TSIZE =>14,
	M6_SCORE =>15,
	M6_GAPS  =>16,  # total gaps length
	M6_POSITIVE =>17,
	M6_QCOV  =>18,
	M6_TCOV  =>19,
};

use constant {
	CT_SYN1122 => 0,
	CT_SYN1122a=> 1,
	CT_1122    => 2, # not synteny
};

my $rerun=0;
my $n;
my $out = $ARGV[0];
my $in_file;
my $out_file;
if (@ARGV<1) {
print <<EOF;
Usage:
  phylogene.pl <PREFIX>
EOF
exit 1;
}
my $clust_dir="$out.cluster";
my $log_dir = "$out.log";
if (! -d $clust_dir) { mkdir $clust_dir; }
if (! -d $log_dir) { mkdir $log_dir; }
my $fin;
my $fout;
my $fout1;
my $fout2;
my $fout3;
my $i=0;
my %species;

use constant {
	SP_CF => "ENSAMX",
	SP_GF => "carAur",
	SP_ZF => "ENSDAR",
	SP_GC => "CTEIDE",
	SP_CC => "CYPCAR",
};

my %copy = (
	CTEIDE => 1,
	CYPCAR => 2,
	ENSAMX => 1,		
	ENSDAR => 1,		
	carAur => 2,		
);

#-----------------------------------------
# write condon.ctl file for PAML codeml
#-----------------------------------------
# seqfile = brown.nuc * sequence data file name
# outfile = mlb * main result file
# treefile = brown.trees * tree structure file name
# noisy = 1 * 0,1,2,3: how much rubbish on the screen
# verbose = 0 * 1: detailed output, 0: concise output
# runmode = 0 * 0: user tree; 1: semi-automatic; 2: automatic
# * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
# * -3: implements a Bayesian method for estimating distance t and 
# *     dN/dS ratio in pairwise comparisons (Angelis et al. 2014). 
# *     The default gamma priors are t ~ G(1.1,1.1) and omega ~ G(1.1, 2.2).
# seqtype = 3 * 1:codons; 2:AAs; 3:codons-->AAs
# CodonFreq = 2 * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
# * ndata = 10
# aaDist = + * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
# * 7:AAClasses
# aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
# * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
# clock = 0 * 0:no clock, 1:clock; 2:local clock
# model = 2
# * models for codons:
# * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
# * models for AAs or codon-translated AAs:
# * 0:poisson, 1:proportional,2:Empirical,3:Empirical+F
# * 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
# NSsites = 0 * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
# * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
# * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
# * 13:3normal>0
# icode = 0 * 0:universal code; 1:mammalian mt; 2-11:see below
# Mgene = 0 * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
# TipDate = 0 100 * TipDate (1) & time unit
# fix_kappa = 0 * 0: estimate kappa; 1: fix kappa at value below
# kappa = 2.5 * initial or fixed kappa
# fix_omega = 0 * 1: omega or omega_1 fixed, 0: estimate
# omega = .4 * initial or fixed omega, for codons or codon-based AAs
# fix_alpha = 1 * 0: estimate gamma shape parameter; 1: fix it at alpha
# alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
# Malpha = 0 * different alphas for genes
# ncatG = 3 * # of categories in dG of NSsites models
# fix_rho = 1 * 0: estimate rho; 1: fix it at rho
# rho = 0. * initial or fixed rho, 0:no correlation
# getSE = 0 * 0: don't want them, 1: want S.E.s of estimates
# RateAncestor = 0 * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
# Small_Diff = .5e-6
# * cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
# * fix_blength = 0 * 0: ignore, -1: random, 1: initial, 2: fixed
# method = 0 * 0: simultaneous; 1: one branch at a time
# * nparK = 0 * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK
# * nhomo = 0 * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2, 5: user
#
# option to set:
# seqfile  outfile  treefile  clock  model  NSite aaDist
my %paml_par = (
seqfile  => "",
outfile  => "",
treefile => "",
noisy    => 1,  # 0,1,2,3: how much rubbish on the screen
verbose  => 0,  # 1: detailed output, 0: concise output
runmode  => 0,  # 0: user tree; 1: semi-automatic; 2: automatic
                # 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
                # -3: implements a Bayesian method for estimating distance t and 
                #     dN/dS ratio in pairwise comparisons (Angelis et al. 2014). 
                #     The default gamma priors are t ~ G(1.1,1.1) and omega ~ G(1.1, 2.2).
seqtype  => 1,   # 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq=> 2,   # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
ndata    => 1,
aaDist   => '0', # 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a, 7:AAClasses
aaRatefile => "wag.dat", # only used for aa seqs with model=empirical(_F)
                         # dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

clock    => 0,  # 0:no clock, 1:clock; 2:local clock
model    => 2,
# models for codons:
# 0:one, 1:b, 2:2 or more dN/dS ratios for branches
# models for AAs or codon-translated AAs:
# 0:poisson, 1:proportional,2:Empirical,3:Empirical+F
# 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
NSsites  => 0, 
# 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
# 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
# 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
# 13:3normal>0
icode    => 0,   # 0:universal code; 1:mammalian mt; 2-11:see below
Mgene    => 0,   # 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
TipDate  => "0 100", # TipDate (1) & time unit
fix_kappa=> 0,   # 0: estimate kappa; 1: fix kappa at value below
kappa    => 2.5, # initial or fixed kappa
fix_omega=> 0,   # 1: omega or omega_1 fixed, 0: estimate
omega    => .4,  # initial or fixed omega, for codons or codon-based AAs
fix_alpha=> 1,   # 0: estimate gamma shape parameter; 1: fix it at alpha
alpha    => 0.,  # initial or fixed alpha, 0:infinity (constant rate)
Malpha   => 0,   # different alphas for genes
ncatG    => 3,   # # of categories in dG of NSsites models
fix_rho  => 1,   # 0: estimate rho; 1: fix it at rho
rho      => 0.,  # initial or fixed rho, 0:no correlation
getSE    => 0,   # 0: don't want them, 1: want S.E.s of estimates
RateAncestor => 0, # (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
Small_Diff => .5e-6,
cleandata   => 0, # remove sites with ambiguity data (1:yes, 0:no)?
fix_blength => 0, # 0: ignore, -1: random, 1: initial, 2: fixed
method      => 0, # 0: simultaneous; 1: one branch at a time
nparK       => undef, # rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK
nhomo       => undef, # 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2, 5: user
);


#-----------------------------------------
# load anchor
#-----------------------------------------
# {{{
$in_file = "$out.anchor";
$fin = base_io::czl_open($in_file, "r") or die "Fail to open $in_file\n";
my %protein;
my %transcript;
my %gene;
while (<$fin>) {
	if (m/^#/) {next;}
	chomp;
	my @t = split "\t",$_,-1;
	my ($sp_id,$sp_name) = split '_', $t[0], 2;
	my ($seq_id,$seq_name) = split '_', $t[1],2;
	my $idx = $t[5];
	my $begin = $t[2];
	my $end   = $t[3];
	my $strand= $t[4];
	my $prot_id = $t[6];
	if ($prot_id =~ /^ENS/) { $prot_id =~ s/\.[0-9]+$//;}
	my $tran_id = $t[7];
	if ($tran_id =~ /^ENS/) { $tran_id =~ s/\.[0-9]+$//;}
	my $gene_id = $t[8];
	if ($gene_id =~ /^ENS/) { $gene_id =~ s/\.[0-9]+$//;}
	my $gene_symbol = $t[9];
	my $to_dup  = $t[10];
	my $pos_idx = $t[11];
	my $clust_id0= $t[12];
	my $clust_size0= $t[13];
	my $clust_id = $t[14];
	my $clust_size = $t[15];
	my $syn_id = $t[16];
	my $nedge  = $t[17];
	if ($nedge==0) { next; }
	$gene{$sp_name}{$gene_id} = {
		idx      => $idx,
		to_dup   => $to_dup,
		seq_name => $seq_name,
		pos_idx  => $pos_idx,
		clust_id0=> $clust_id0,
		clust_size0 => $clust_size0,
		clust_id    => $clust_id,
		clust_size  => $clust_size,
		syn_id      => $syn_id,
		protein_id    => $prot_id,
		transcript_id => $tran_id,
		gene_id       => $gene_id,
		name          => $gene_symbol,
		het      => 0,
	};
	if ($prot_id ne "" && $prot_id ne ".") {
		$protein{$sp_name}{$prot_id}{gene_id} = $gene_id;
		$protein{$sp_name}{$prot_id}{transcript_id} = $tran_id;
	}
	if ($tran_id ne "" && $tran_id ne ".") {
		$transcript{$sp_name}{$tran_id}{gene_id}    = $gene_id;
		$transcript{$sp_name}{$tran_id}{protein_id} = $prot_id;
	}
	$species{$sp_name}{n_gene}++;
}
close $fin;
my @sps = sort keys(%species);
foreach my $sp (keys(%gene)) {
	foreach my $id (keys(%{$gene{$sp}})) {
		my $g = $gene{$sp}{$id};
		if (!defined $g->{to_dup} || $g->{to_dup} eq "-1" || $g->{to_dup} eq ".") {  }
		else {
			$g->{to_dup} = $protein{$sp}{$g->{to_dup}}{gene_id};
			push @{$gene{$sp}{$g->{to_dup}}{dup}}, [ $id, $g->{transcript_id}, $g->{protein_id}, $g->{name} ];
#			delete $gene{$sp}{$id};
		}
	}
}
foreach my $sp (keys(%gene)) {
	foreach my $id (keys(%{$gene{$sp}})) {
		my $g = $gene{$sp}{$id};
		if (!exists $g->{dup}) { $g->{dup}=[]; }
	}
}
# }}}

# load anchor annotation
# $in_file = "/data/chenz11/goldfish/11549472/sergey_canu70x/arrow/carAur03/ens85_GC_CC_GF.bed";
#$in_file = "all5.6.bed";
#$fin = base_io::czl_open($in_file, "r") or die "Fail to open $in_file\n";
#while (<$fin>) {
#	if (m/^#/) {next;}
#	chomp;
#	my @t = split "\t",$_,-1;
#	my ($sp, $gene_id,$transcript_id,$protein_id,$name,$transcript_len,$protein_len) = split /\|/, $t[3];
#	if (defined $protein_id && $protein_id ne ".") {
#		if (!exists $protein{$sp}{$protein_id}) { next; }
#		$transcript{$sp}{$transcript_id}{gene_id} = $gene_id;
#		$transcript{$sp}{$transcript_id}{begin}   = $t[1];
#		$transcript{$sp}{$transcript_id}{end}     = $t[2];
#		$transcript{$sp}{$transcript_id}{strand}  = $t[5];
#		$transcript{$sp}{$transcript_id}{name}    = $name;
#		$transcript{$sp}{$transcript_id}{size}    = $transcript_len;
#		$transcript{$sp}{$transcript_id}{protein_id} = $protein_id;
#		$protein{$sp}{$protein_id}{transcript_id} = $transcript_id;
#		$protein{$sp}{$protein_id}{gene_id} = $gene_id;
#		$protein{$sp}{$protein_id}{name}    = $name;
#		$protein{$sp}{$protein_id}{begin}   = $t[1];
#		$protein{$sp}{$protein_id}{end}     = $t[2];
#		$protein{$sp}{$protein_id}{strand}  = $t[5];
#	}
#}
#close $fin;

# -----------------------------------------------
# set cluster names and size
# -----------------------------------------------
# {{{
my @clust;
foreach my $sp (keys(%gene)) {
	foreach my $id (keys(%{$gene{$sp}})) {
		my $aa = $gene{$sp}{$id};
		my $tid = $aa->{transcript_id};
		my $cid = $aa->{clust_id};
		if ($cid>=0) {
			push @{$clust[$cid]{anchor}{$sp}}, $id;
			$clust[$cid]{size}++;
		}
		$aa->{sclust_id} = -1;
	}
}
for (my $cid=0; $cid<@clust; $cid++) {
	my $c = $clust[$cid];
	my $sp = "ENSDAR";
	my $name;
	my @name;
	if (!exists $c->{anchor}{$sp}) { $c->{anchor}{$sp}=[]; }
	if (@{$c->{anchor}{$sp}}>0) {
		my $bb = $c->{anchor}{$sp};
		my $id = $bb->[0];
		if (exists $gene{$sp}{$id}) {
			my $p = $gene{$sp}{$id};
			if (defined $p->{name} && $p->{name} ne ".") {push @name, $p->{name};}
		}
	}
	$c->{names}=\@name;
	if (@name == 0) {
		$c->{name} = ".";
	} else {
		$c->{name} = join ",", @name;
	}
}
$out_file = "$out.cluster_size.txt";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
print $fout "#ClusterID\tGeneName\tSize", join("\t",@sps), "\n";
for (my $cid=0; $cid<@clust; $cid++) {
	my $c = $clust[$cid];
	if (!defined $c) { next; }
	my $aa = $c->{anchor};
	print $fout "$cid\t$c->{name}\t$c->{size}";
	foreach my $sp (@sps) {
		if (exists $aa->{$sp}) { print $fout "\t", $#{$aa->{$sp}}+1; }
		else { print $fout "\t0"; }
	}
	print $fout "\n";
}
close $fout;
# }}}

# cat out11.cluster_size.txt | tail -n +2 | awk '$3>=4 {print $7"__"$8"__"$5}' | sort | uniq -c | awk '{print $2"\t"$1}' | sed 's/__/\t/g' > out11.cluster_size.hist.txt

#------------------------------------------
# load contig depth
#------------------------------------------
$in_file = "carAur01.ctgInfo2";
my %contig;
$fin = base_io::czl_open($in_file, "r") or die "Fail to open $in_file\n";
while (<$fin>) {
	chomp;
	my @t=split "\t",$_,-1;
	$contig{$t[1]}{rdepth} = $t[9];
}
close $fin;
# set contig depth for anchor
{
	my $sp = "carAur";
	foreach my $id (keys(%{$gene{$sp}})) {
		my $aa = $gene{$sp}{$id};
		my $chr = $aa->{seq_name};
		if (exists $contig{$chr}{rdepth}) {
			$aa->{rdepth} = $contig{$chr}{rdepth};
		}
	}
}

# -----------------------------------------------
# count homologs for each chromosome pairs (zf,cc) or (zf,gf)
# -----------------------------------------------
# {{{
my %chr_pair;
my %chr;
my %chr_nums;
{
	foreach (my $cid=0; $cid<@clust; $cid++) {
		my $c= $clust[$cid];
		if (!defined $c || $c->{size}<=2) { next; }
		if (!exists $c->{anchor}{ENSDAR}) { next; }
		if (@{$c->{anchor}{ENSDAR}}>5) { next; }
		my @cc_cid;
		foreach my $zf_id (@{$c->{anchor}{ENSDAR}}) {
			my $zf_seq_name = $gene{ENSDAR}{$zf_id}{seq_name};
			my $chr = $zf_seq_name;
			$chr{ENSDAR}{$chr}++;
			$chr=~ s/^chr//;
			if ($chr !~ /^[0-9]/) { next; }
			foreach my $sp (qw(CTEIDE CYPCAR carAur)) {
				if (exists $c->{anchor}{$sp} && @{$c->{anchor}{$sp}}<=5*$copy{$sp}) {
					foreach my $id2 (@{$c->{anchor}{$sp}}) {
						my $g    = $gene{$sp}{$id2};
						my $chr2 = $g->{seq_name};
						if ($chr2=~m/LG([0-9]+)/) {
							$chr2=$1;
							$chr_nums{ENSDAR}{$chr}++;
							$chr_nums{$sp}{$chr2}++;
							$chr_pair{ENSDAR}{$sp}{$chr}{$chr2}++;
							$chr_pair{$sp}{ENSDAR}{$chr2}{$chr}++;
							$chr{$sp}{$chr2}++;
						} elsif ($chr2=~m/chr([0-9]+)/) {
							$chr2=$1;
							$chr_nums{ENSDAR}{$chr}++;
							$chr_nums{$sp}{$chr2}++;
							$chr_pair{ENSDAR}{$sp}{$chr}{$chr2}++;
							$chr_pair{$sp}{ENSDAR}{$chr2}{$chr}++;
							$chr{$sp}{$chr2}++;
						} else {
							$chr_pair{ENSDAR}{$sp}{$chr}{$chr2}++;
							$chr_pair{$sp}{ENSDAR}{$chr2}{$chr}++;
						}
					}
				}
			}
		}

		if (!exists $c->{anchor}{CYPCAR} || @{$c->{anchor}{CYPCAR}}>5*$copy{CYPCAR}) { next; }
		if (!exists $c->{anchor}{carAur} || @{$c->{anchor}{carAur}}>5*$copy{carAur}) { next; }
		foreach my $cc_id (@{$c->{anchor}{CYPCAR}}) {
			my $cc_seq_name = $gene{CYPCAR}{$cc_id}{seq_name};
			my $chr = $cc_seq_name;
			if ($chr=~ s/^LG([0-9]+)//) { $chr = $1; }
			else { next; }
			foreach my $gf_id (@{$c->{anchor}{carAur}}) {
				my $gf_p   = $gene{carAur}{$gf_id};
				my $gf_chr = $gf_p->{seq_name};
				if ($gf_chr=~/LG([0-9]+)/) {
					$gf_chr=$1;
					$chr_pair{CYPCAR}{carAur}{$chr}{$gf_chr}++;
				}
			}
		}
	}
#	foreach my $sp (qw(CTEIDE CYPCAR carAur)) {
#		foreach my $chr (keys($chr_pair{$sp}{ENSDAR})) {
#		}
#	}

	foreach my $sp (qw(CTEIDE CYPCAR carAur)) {
		my @chrs = sort {$a<=>$b} keys %{$chr_nums{ENSDAR}};
		my @chrs_out = @chrs;
		foreach my $chr (@chrs_out) { $chr = "LG$chr"; }
		$out_file = "$out.ENSDAR_$sp.chr_homolog_count.txt";
		$fout1 = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
		print $fout1 join("\t", @chrs_out), "\n";
		foreach my $chr1 (@chrs) {
			my $aa = $chr_pair{ENSDAR}{$sp}{$chr1};
			print $fout1 "chr$chr1";
			foreach my $chr2 (@chrs) {
				if (!exists $aa->{$chr2}) { $aa->{$chr2}=0; }
				print $fout1 "\t", $aa->{$chr2};
			}
			print $fout1 "\n";
		}
		close $fout1;
	}

	{
		$out_file = "$out.CYPCAR_carAur.chr_homolog_count.txt";
		$fout1 = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
		my @chrs = sort {$a<=>$b} keys(%{$chr{carAur}});
		my @chrs_out = @chrs;
		print $fout1 join("\t", @chrs_out), "\n";
		foreach my $chr (@chrs_out) { $chr = "LG$chr"; }
		foreach my $chr1 (sort {$a<=>$b} keys(%{$chr_pair{CYPCAR}{carAur}})) {
			print $fout1 "LG$chr1";
			my $aa=$chr_pair{CYPCAR}{carAur}{$chr1};
			foreach my $chr2 (@chrs) {
				if (!exists $aa->{$chr2}) { $aa->{$chr2}=0; }
				print $fout1 "\t", $aa->{$chr2};
			}
			print $fout1 "\n";
		}
	}
	close $fout1;
}

my %match_chr;

foreach my $sp (qw(CTEIDE CYPCAR carAur)) {
	foreach my $chr (keys(%{$chr_pair{ENSDAR}{$sp}})) {
		my $aa = $chr_pair{ENSDAR}{$sp}{$chr};
		my @chrs = sort {$aa->{$b} <=> $aa->{$a} } keys(%{$aa});
		my $n = $aa->{$chrs[0]};
		my $i;
		for ($i=0; $i<@chrs && $i<$copy{$sp}; $i++) {
			push @{$match_chr{$sp}{$chr}}, $chrs[$i];
		}
		if ($i<@chrs) {
			foreach my $chr2 (@chrs[$i..$#chrs]) {
				if ($aa->{$chr2}>100 || $aa->{$chr2}>$n*0.8) {
					push @{$match_chr{$sp}{$chr}}, $chr2;
				} else {
					last;
				}
			}
		}
	}
	$out_file = "$out.chr_pairs.ENSDAR.$sp.txt";
	$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
	foreach my $chr (sort {$a <=> $b} keys(%{$match_chr{$sp}})) {
		print $fout "$chr\t", join(",", @{$match_chr{$sp}{$chr}}), "\n";
	}
	close $fout;
}
# }}}

# -----------------------------------------------
# load edges
# -----------------------------------------------
# {{{
my %link;
if (-f "all5.f3.join.f.edge") {
	$in_file = "all5.f3.join.f.edge";
	$fin = base_io::czl_open($in_file, "r") or die "Fail to open $in_file\n";
	while (<$fin>) {
		if (m/^#/) { next; }
		chomp;
		my @t = split "\t",$_,-1;
		my ($sp1,$id1) = split /\|/, $t[0];
		my ($sp2,$id2) = split /\|/, $t[1];
		my $score = $t[3];
		my $iden  = $t[4];
		if (exists $link{$sp2}{$sp1}{$id2}{$id1}) {
			if ($link{$sp2}{$sp1}{$id2}{$id1}[1] < $iden) {
				$link{$sp1}{$sp2}{$id1}{$id2}= [ $score, $iden ];
				$link{$sp2}{$sp1}{$id2}{$id1}= [ $score, $iden ];
			}
		} else {
			$link{$sp1}{$sp2}{$id1}{$id2}= [ $score, $iden ];
			$link{$sp2}{$sp1}{$id2}{$id1}= [ $score, $iden ];
		}
	}
	close $fin;
}
else {
	$in_file = "all5.f3.join.edge";
	$fin = base_io::czl_open($in_file, "r") or die "Fail to open $in_file\n";
	my $rm=0;
	while (<$fin>) {
		if (m/^#/) { next; }
		chomp;
		my @t = split "\t",$_,-1;
		my ($sp1,$id1) = split /\|/, $t[0];
		my ($sp2,$id2) = split /\|/, $t[1];
		if (!exists $protein{$sp1}{$id1}) { $rm++; next; }
		if (!exists $protein{$sp2}{$id2}) { $rm++; next; }
		$id1=$protein{$sp1}{$id1}{gene_id};
		$id2=$protein{$sp2}{$id2}{gene_id};
		if ($gene{$sp1}{$id1}{to_dup} ne "." && $gene{$sp1}{$id1}{to_dup}!~/^\s$/) {
			$id1 = $gene{$sp1}{$id1}{to_dup};
		}
		if ($gene{$sp2}{$id2}{to_dup} ne "." && $gene{$sp2}{$id2}{to_dup}!~/^\s$/) {
			$id2 = $gene{$sp2}{$id2}{to_dup};
		}
		my $score = $t[3];
		my $iden  = $t[4];
		if (exists $link{$sp2}{$sp1}{$id2}{$id1}) {
			if ($link{$sp2}{$sp1}{$id2}{$id1}[1] < $iden) {
				$link{$sp1}{$sp2}{$id1}{$id2}= [ $score, $iden ];
				$link{$sp2}{$sp1}{$id2}{$id1}= [ $score, $iden ];
			}
		} else {
			$link{$sp1}{$sp2}{$id1}{$id2}= [ $score, $iden ];
			$link{$sp2}{$sp1}{$id2}{$id1}= [ $score, $iden ];
		}
	}
	close $fin;
	print "Remove Edge Num = $rm\n";
}

my %link_iden_max;
my %link_best;
print "#---------------------------\n";
print "# Good BLAST link degree counts\n";
foreach my $sp1 (keys(%link)) {
my $n1 = $copy{$sp1};
foreach my $sp2 (keys(%{$link{$sp1}})) {
my $n2 = $copy{$sp2};
if ($n1>$n2 && $n2!=1) { $n1/=$n2; $n2=1; }
elsif ($n1<$n2 && $n1!=1) { $n2/=$n1; $n1=1; }
my @count = (0) x 11;
foreach my $id1 (keys(%{$link{$sp1}{$sp2}})) {
	my $aa = $link{$sp1}{$sp2}{$id1};
	my @id2s = sort {$aa->{$b}[1] <=> $aa->{$a}[1]} keys(%$aa);
	my %id2h;
	my $iden_max = $aa->{$id2s[0]}[1];
	if ($sp1 eq $sp2 && $id1 eq $id2s[0]) {
		if (@id2s>=2) { $iden_max=$aa->{$id2s[1]}[1]; }
		else { $iden_max=0; }
	}
	$link_iden_max{$sp1}{$sp2}{$id1} = $iden_max;
	for (my $i=0; $i<@id2s && $i<=5*$copy{$sp2}; $i++) {
		my $id2 = $id2s[$i];
		if ($sp1 eq $sp2 && $id1 eq $id2) { next; }
		if ( $aa->{$id2}[1]>=0.9*$iden_max || @{$link_best{$sp1}{$sp2}{$id1}}<$copy{$sp2}) {
			push @{$link_best{$sp1}{$sp2}{$id1}}, $id2;
			$id2h{$id2}++;
		}
		elsif ( $aa->{$id2}[1]>=0.75*$iden_max ) { $id2h{$id2}++; }
		else { last; }
	}
	foreach my $id2 (keys(%$aa)) {
		if ($sp1 eq $sp2 && $id1 eq $id2) { next; }
		if (!exists $id2h{$id2}) { delete $aa->{$id2}; }
	}
	$count[scalar keys(%$aa)]++;
}
print "$sp1\t$sp2\t", join("\t",@count), "\n";
} }
print "#---------------------------\n";

if (! -f "all5.f3.join.f.edge") {
	$out_file = "all5.f3.join.f.edge";
	$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
	foreach my $sp1 (sort keys(%link)) {
	foreach my $sp2 (sort keys(%{$link{$sp1}})) {
	if ($sp1 gt $sp2) { next; }
	foreach my $id1 (sort keys(%{$link{$sp1}{$sp2}})) {
	my $aa = $link{$sp1}{$sp2}{$id1};
	my @id2s = sort {$aa->{$b}[1] <=> $aa->{$a}[1]} keys(%$aa);
	foreach my $id2 (@id2s) {
		my $bb = $link{$sp1}{$sp2}{$id1}{$id2};
		print $fout "$sp1|$id1\t$sp2|$id2\t+\t$bb->[0]\t$bb->[1]\n";
	} } } }
	close $fout;
}
# }}}
#-----------------------------------
foreach my $sp (keys(%gene)) {
	foreach my $id (keys(%{$gene{$sp}})) {
		my $g = $gene{$sp}{$id};
		if ($g->{to_dup} ne "" && $g->{to_dup} ne ".") { delete $gene{$sp}{$id}; }
	}
}



#my %zf_to;
#for (my $i=1; $i<=25; $i++) { $zf_to{CYPCAR}[$i] = [$i*2-1, $i*2]; }
#for (my $i=1; $i<=25; $i++) { $zf_to{carAur}[$i] = [$i, $i+25]; }
#$zf_to{CTEIDE} = [ [], [8], [22], [2,13], [11], [17], [18], [1], [21], [16], [24], [10], [9], [5], [13], [15], [12,19], [23], [7], [4], [14], [3], [24], [6], [19], [20] ];
# -----------------------------------------------
# fetch ortholog-ohnolog cluster
# -----------------------------------------------
my @sclust;
my $n_het=0;
{
	my @gf_n= (0) x 4;
	my @cc_n= (0) x 4;
	my $n=0;
	my $scid=0;
# add anchor on LG
# {{{
	foreach (my $cid=0; $cid<@clust; $cid++) {
		my $c= $clust[$cid];
		if (!defined $c) { next; }
		if ($c->{size}<=2) { next; }
		if (!exists $c->{anchor}{ENSDAR}) { next; }
		my %sz;
		my $large=0;
		foreach my $sp (@sps) {
			if (exists $c->{anchor}{$sp}) {$sz{$sp} = scalar @{$c->{anchor}{$sp}};}
			else { $sz{$sp} = 0; }
			if ($sp eq 'CYPCAR' || $sp eq 'carAur') {
				if ($sz{$sp}>10) {$large++;}
			} else {
				if ($sz{$sp}>5) {$large++;}
			}
		}
		if ($large) { next; };
		my $n_zf= $sz{ENSDAR};
		my $n_cc= $sz{CYPCAR};
		my $n_gf= $sz{carAur};
		my $n_gc= $sz{CTEIDE};
		foreach my $zf_id (@{$c->{anchor}{ENSDAR}}) {
			my $chr = $gene{ENSDAR}{$zf_id}{seq_name};
			$chr=~ s/^chr//;
			if ($chr !~ /^[0-9]/) { next; }
			my @gf_seq_name = ( "LG$chr", "LG" . ($chr+25) );
			my @cc_seq_name = ( "LG".($chr*2-1), "LG".($chr*2) );
			my $gf_ids;
			my $cc_ids;

			my %sc;
			foreach my $sp (@sps) { $sc{anchor}{$sp} = []; }
			my %lg;
			$lg{ENSDAR}{$chr} = [$zf_id];
			push @{$sc{anchor}{ENSDAR}}, $zf_id;
			$sc{size} = 1;
			foreach my $sp2 (qw(CTEIDE CYPCAR carAur)) {
				$sc{anchor}{$sp2}=[];
				foreach my $chr2 (@{$match_chr{$sp2}{$chr}}) {
#					if ($chr2=~m/^[0-9XYZW][0-9]*$/) {
						$lg{$sp2}{$chr2} = [];
#					}
				}
				if (exists $c->{anchor}{$sp2}) {
					foreach my $id2 (@{$c->{anchor}{$sp2}}) {
						my $chr2 = $gene{$sp2}{$id2}{seq_name};
						if ($chr2=~m/LG([0-9]+)/)  {$chr2=$1; }
						if (exists $lg{$sp2}{$chr2}) {
							my $bb = $lg{$sp2}{$chr2};
							if ( @$bb==0 ) {
								push @{$sc{anchor}{$sp2}}, $id2;
								$sc{size}++;
							}
							push @{$lg{$sp2}{$chr2}}, $id2;
						}
					}
				}
			}

			if ($sc{size}<3) { next; }

			if ( @{$sc{anchor}{CYPCAR}}>=2 ) {
				$cc_n[0]++;
			} elsif ( @{$sc{anchor}{CYPCAR}}==1 ) {
				$cc_n[1]++;
			}
			if ( @{$sc{anchor}{carAur}}>=2 ) {
				$gf_n[0]++;
			} elsif ( @{$sc{anchor}{carAur}}==1 ) {
				$gf_n[1]++;
			}
			if ( @{$sc{anchor}{CYPCAR}}>=2 && @{$sc{anchor}{carAur}}>=2 ) { $n++; }

			$sc{name} = $gene{ENSDAR}{$zf_id}{name};
			$sc{clust_id} = $cid;
			foreach my $sp (keys(%{$sc{anchor}})) {
				foreach my $id (@{$sc{anchor}{$sp}}) {
					$gene{$sp}{$id}{sclust_id} = $scid;
				}
			}
			push @sclust, \%sc;
			push @{$c->{sclust_ids}}, $scid;
			$scid++;
		}
	}
	print "N chr-synteny cluster = " . ($#sclust+1) . "\n";
	print "N_chr_syn_gf = ($gf_n[0],$gf_n[1])\n";
	print "N_chr_syn_cc = ($cc_n[0],$cc_n[1])\n";
	print "N_chr_syn_gf_cc = $n\n";
# }}}

# add anchor not on LG but in same cluster
# {{{
	$n=0;
	my $nc2=0;
	for (my $scid=0; $scid<@sclust; $scid++) {
		my $sc = $sclust[$scid];
		my $cid = $sc->{clust_id};
		my $c  = $clust[$cid];
		my $sca = $sc->{anchor};
		my $ca  = $c->{anchor};
		my $zf_id = $sca->{ENSDAR}[0];
		my $zf_chr = $gene{ENSDAR}{$zf_id}{seq_name};
		my %n0;
		my %n1;
		my %anchor;
		foreach my $sp (@sps) {
			$n0{$sp}=0;
			$n1{$sp}=0;
			$anchor{$sp} = [];
			foreach my $id (@{$ca->{$sp}}) {
				my $g = $gene{$sp}{$id};
				if (exists $g->{het} && $g->{het}!=0) { next; }
				if ($g->{sclust_id}>=0) { $n0{$sp}++; }
				else {
					push @{$anchor{$sp}}, $id;
					$n1{$sp}++;
				}
			}
		}
		if ($n0{ENSDAR}>1) { $nc2++; }
		foreach my $sp (qw(ENSAMX CTEIDE CYPCAR carAur)) {
			my $m0 = @{$sca->{$sp}};
			my $m1 = $copy{$sp} - $m0;
			if ($m0<$copy{$sp}) {
				if ($n1{$sp}>0) {
					if ($n1{$sp} + $m0 <= $copy{$sp}) {
						foreach my $id (@{$anchor{$sp}}) { $gene{$sp}{$id}{sclust_id} = $scid; }
						push @{$sca->{$sp}}, @{$anchor{$sp}};
					} else {
						my $aa = $link{ENSDAR}{$sp}{$zf_id};
						@{$anchor{$sp}} = sort { 
							if (exists $aa->{$a}) {
								if (exists $aa->{$b}) {
									-($aa->{$b}[1] <=> $aa->{$a}[1]);
								} else {
									return -1;
								}
							} else {
								if (exists $aa->{$b}) {
									return 1;
								} else {
									return 0;
								}
							}
						} @{$anchor{$sp}};
						if ($sp eq 'carAur') {
							for (my $i=0; $i<@{$anchor{$sp}}; $i++) {
								my $id1 = $anchor{$sp}[$i];
								my $bad=0;
								my $g1 = $gene{$sp}{$id1};
								my $chr1 = $g1->{seq_name};
								if (defined $g1->{rdepth} && $g1->{rdepth}<0.8) {
									foreach my $id2 (@{$sca->{$sp}}) {
										if (exists $link{$sp}{$sp}{$id1}{$id2} && $link{$sp}{$sp}{$id1}{$id2}[1]>98) { $bad=1; $g1->{het}=1; $n_het++;}
									}
								}
								if ($bad) { next; }
								push @{$sca->{$sp}}, $id1;
								$gene{$sp}{$id1}{sclust_id} = $scid;
								if (@{$sca->{$sp}}>=$copy{$sp}) { last; }
							}
						} else {
							for (my $i=0; $i<$m1; $i++) {
								my $id = $anchor{$sp}[$i];
								push @{$sca->{$sp}}, $id;
								$gene{$sp}{$id}{sclust_id} = $scid;
							}
						}
					}
				}
			}
		}

		if ( @{$sca->{CYPCAR}}>=2 ) {
			$cc_n[2]++;
		} elsif ( @{$sca->{CYPCAR}}==1 ) {
			$cc_n[3]++;
		}
		if ( @{$sca->{carAur}}>=2 ) {
			$gf_n[2]++;
		} elsif ( @{$sca->{carAur}}==1 ) {
			$gf_n[3]++;
		}
		if ( @{$sca->{CYPCAR}}>=2 && @{$sca->{carAur}}>=2 ) { $n++; }
	}
	print "N chr-synteny cluster = " . ($#sclust+1) . "\n";
	print "N_chr_syn_gf = ($gf_n[2],$gf_n[3])\n";
	print "N_chr_syn_cc = ($cc_n[2],$cc_n[3])\n";
	print "N_chr_syn_gf_cc = $n\n";
# }}}

	$n=0;
	$cc_n[2]=$cc_n[3]=0;
	$gf_n[2]=$gf_n[3]=0;
	for (my $scid=0; $scid<@sclust; $scid++) {
		my $sc = $sclust[$scid];
		my $cid = $sc->{clust_id};
		my $c  = $clust[$cid];
		my $sca = $sc->{anchor};
		my $ca  = $c->{anchor};
		my $zf_id = $sca->{ENSDAR}[0];
		my $zf_chr = $gene{ENSDAR}{$zf_id}{seq_name};
		foreach my $sp (qw(ENSAMX CTEIDE CYPCAR carAur)) {
			my $m0 = @{$sca->{$sp}};
			my $m1 = $copy{$sp} - $m0;
			if ($m0<$copy{$sp}) {
				my $aa = $link{ENSDAR}{$sp}{$zf_id};
				my @ids = keys(%$aa);
				@ids = sort {$aa->{$b}[1]<=>$aa->{$a}[1]} @ids; # sort by identity
				foreach my $id1 (@ids) {
					my $g1 = $gene{$sp}{$id1};
					if ($g1->{sclust_id}>=0) { next; }
					if ($g1->{het}!=0) { next; }
					my $cid1 = $g1->{clust_id};
					my $c1 = $clust[$cid1];
					my $ca1 = $clust[$cid1]{anchor};
					my $chr1 = $g1->{seq_name};
					if ($chr1=~m/(chr|LG)([0-9]+)/) { $chr1 = $2; }
					if (@{$ca1->{ENSDAR}}==0) {
						if ($sp eq 'carAur') {
							my $bad=0;
							if (defined $g1->{rdepth} && $g1->{rdepth}<0.8) {
								foreach my $id2 (@{$sca->{$sp}}) {
									my $g2 = $gene{$sp}{$id2};
									my $chr2 = $g2->{seq_name};
									if ($chr2=~m/(chr|LG)([0-9]+)/) { $chr2 = $2; }
#									if ($chr1 eq $chr2) { next; }
									if (exists $link{$sp}{$sp}{$id1}{$id2} && $link{$sp}{$sp}{$id1}{$id2}[1]>98) { $bad=1; $g1->{het}=1; $n_het++;}
								}
							}
							if ($bad) { next; }
							push @{$sca->{$sp}}, $id1;
							$gene{$sp}{$id1}{sclust_id} = $scid;
							if (@{$sca->{$sp}}>=$copy{$sp}) { last; }
						} else {
							push @{$sca->{$sp}}, $id1;
							$gene{$sp}{$id1}{sclust_id} = $scid;
						}
						if (@{$sca->{$sp}}==$copy{$sp}) { last; }
					}
				}
			}
		}
		if ( @{$sca->{CYPCAR}}>=2 ) {
			$cc_n[2]++;
		} elsif ( @{$sca->{CYPCAR}}==1 ) {
			$cc_n[3]++;
		}
		if ( @{$sca->{carAur}}>=2 ) {
			$gf_n[2]++;
		} elsif ( @{$sca->{carAur}}==1 ) {
			$gf_n[3]++;
		}
		if ( @{$sca->{CYPCAR}}>=2 && @{$sca->{carAur}}>=2 ) { $n++; }
	}
	print "N chr-synteny cluster = " . ($#sclust+1) . "\n";
	print "N_chr_syn_gf = ($gf_n[2],$gf_n[3])\n";
	print "N_chr_syn_cc = ($cc_n[2],$cc_n[3])\n";
	print "N_chr_syn_gf_cc = $n\n";

	print "N_heterozygo_gf = $n_het\n";
}

#-----------------------------------------------
# filter goldfish heterozygous TODO
#-----------------------------------------------

#------------------------------
# Output scluster size
#------------------------------
$out_file = "$out.scluster_size.txt";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
print $fout "#ClusterID\tGeneName\tSize\t", join("\t",@sps), "\n";
for (my $scid=0; $scid<@sclust; $scid++) {
	my $sc = $sclust[$scid];
	my $sca = $sc->{anchor};
	print $fout "$scid\t$sc->{name}\t$sc->{size}\t$sc->{clust_id}";
	foreach my $sp (@sps) {
		print $fout "\t", $#{$sca->{$sp}}+1;
		print $fout "\t";
		my $j=0;
		if (@{$sca->{$sp}}==0) { print "."; }
		else {
			foreach my $id (@{$sca->{$sp}}) {
				if ($j>0) { print $fout ","; }
				my $g = $gene{$sp}{$id};
				print $fout "$g->{seq_name}:$id";
			}
		}
	}
	print $fout "\n";
}
close $fout;

{
#-----------------------------------------------
# WGD gene lost TODO
#-----------------------------------------------
	my %chr_size;
	my %chr_size_scid;
	for (my $scid=0; $scid<@sclust; $scid++) {
		my $sc = $sclust[$scid];
		my $sca = $sc->{anchor};
		my $zf_id = $sca->{ENSDAR}[0];
		my $zf_g  = $gene{ENSDAR}{$zf_id};
		my $zf_chr= $zf_g->{seq_name};
		if ($zf_chr =~ /^chr([0-9]+)/) { $zf_chr=$1; }
		my $cc_sz = @{$sca->{CYPCAR}};
		my $gf_sz = @{$sca->{carAur}};
		my $k=$cc_sz*3+$gf_sz;
		my $key;
		if ($gf_sz==2) {
			my $gf_id0 = $sca->{carAur}[0];
			my $gf_g0  = $gene{carAur}{$gf_id0};
			my $gf_chr0= $gf_g0->{seq_name};
			my $gf_id1 = $sca->{carAur}[1];
			my $gf_g1  = $gene{carAur}{$gf_id1};
			my $gf_chr1= $gf_g1->{seq_name};
			if ($gf_chr0 =~ m/(chr|LG)[0-9]+/) {
				if ($gf_chr1 =~ m/(chr|LG)[0-9]+/) {
					$key = join(":", sort ($gf_chr0,$gf_chr1));
				} else {
					$key = "$gf_chr0:CTG";
				}
			} else {
				if ($gf_chr1 =~ m/(chr|LG)[0-9]+/) {
					$key = "$gf_chr1:CTG";
				} else {
					$key = "CTG:CTG";
				}
			}
		} elsif ($gf_sz==1) {
			my $gf_id = $sca->{carAur}[0];
			my $gf_g  = $gene{carAur}{$gf_id};
			my $gf_chr= $gf_g->{seq_name};
			if ($gf_id =~ m/(chr|LG)[0-9]+/) {
				$key = "$gf_chr:.";
			} else {
				$key = "CTG:.";
			}
		} else {
			$key = ".:.";
		}
		if (!exists $chr_size{$zf_chr}{$key}) {
			$chr_size{$zf_chr}{$key} = [ (0) x 9 ];
			for (my $k=0; $k<9; $k++) {
				$chr_size_scid{$zf_chr}{$key}[$k] = [];
			}
		}
		$chr_size{$zf_chr}{$key}[$k]++;
		push @{$chr_size_scid{$zf_chr}{$key}[$k]}, $scid;
	}

	$out_file = "$out.scluster_size_zf_chr_count.txt";
	$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
	print $fout "ENSDAR.CHROM\tcarAur.CHROM";
	for (my $i=0; $i<3; $i++) {
		for (my $j=0; $j<3; $j++) {
			print $fout "\t$i:$j";
		}
	}
	for (my $i=0; $i<3; $i++) {
		for (my $j=0; $j<3; $j++) {
			print $fout "\t$i:$j.sclusters";
		}
	}
	print $fout "\n";
	foreach my $chr (sort {$a<=>$b} keys(%chr_size)) {
	foreach my $key2 (keys(%{$chr_size{$chr}})) {
		print $fout "$chr\t$key2";
		my $aa = $chr_size{$chr}{$key2};
		for (my $i=0; $i<3; $i++) {
			for (my $j=0; $j<3; $j++) {
				my $k = $i*3+$j;
				print $fout "\t$aa->[$k]";
			}
		}
		for (my $i=0; $i<3; $i++) {
			for (my $j=0; $j<3; $j++) {
				print $fout "\t";
				my $k = $i*3+$j;
				my $aa = $chr_size_scid{$chr}{$key2}[$k];
				if (@$aa==0) {
					print $fout ".";
				} else {
					for (my $i1=0; $i1<@$aa; $i1++) {
						my $scid1 = $aa->[$i1];
						my $sc1 = $sclust[$scid1];
						if ($i1>0) { print $fout ","; }
						print $fout "$scid1:$sc1->{name}";
					}
				}
			}
		}
		print $fout "\n";
	}
	}
	close $fout;
}


#----------------------------------------------


# --DELETE-- 20171128
##------------------------------------------
## load similarity from blast M6 format
##------------------------------------------
#	$in_file = "carAur.carAur.f2.join.m6.gz";
#	$fin = base_io::czl_open($in_file, "r") or die "Fail to open $in_file\n";
#	my %align;
#	while (<$fin>) {
#		chomp;
#		my @t=split "\t",$_,-1;
#		my $piden = $t[M6_PIDEN];
#		my $tid0 = $t[M6_TID];
#		my $qid0 = $t[M6_QID];
#		my ($tsp,$tid) = split /\|/, $tid0;
#		my ($qsp,$qid) = split /\|/, $qid0;
#		if ($tsp eq "carAur" && $qsp eq "carAur" && $tid0 ne $qid0) {
#			$align{$tsp}{$qsp}{$tid}{$qid} = \@t;
#		}
#	}
#	close $fin;
## ------------------------------------------
#
#for (my $cid=0; $cid<@clust; $cid++) {
#	my $c = $clust[$cid];
#	my $sp = "carAur";
#	my $name = ".";
#	if (exists $c->{anchor}{$sp} && @{$c->{anchor}{$sp}}>1) {
#		my $bb = $c->{anchor}{$sp};
#		for (my $i=0; $i<@$bb-1; $i++) {
#			my $pid1 = $bb->[$i][0];
#			my $chr1 = $protein{$sp}{$pid1}{seq_name};
#			my $p1 = $protein{$sp}{$pid1};
#			for (my $j=$i+1; $j<@$bb; $j++) {
#				my $pid2 = $bb->[$j][0];
#				my $chr2 = $protein{$sp}{$pid2}{seq_name};
#				my $p2 = $protein{$sp}{$pid2};
#				my $aa;
#				if (exists $align{$sp}{$sp}{$pid1}{$pid2}) {$aa=$align{$sp}{$sp}{$pid1}{$pid2};}
#				elsif (exists $align{$sp}{$sp}{$pid2}{$pid1}) {$aa=$align{$sp}{$sp}{$pid2}{$pid1};}
#				if ($chr1 =~ m/^LG/ && $chr2 =~ m/^LG/) { next; }
#				if ($chr1 !~ m/^LG/ && $contig{$chr1}{rdepth}<0.75 && !exists $p2->{het}) {
#					if (!defined $aa) {
#						$p1->{het} = 1;
#					} else {
#						my $piden = $aa->[M6_PIDEN];
#						if ($piden>97) {$p1->{het} = 2;}
#					}
#				}
#				if ($chr2 !~ m/^LG/ && $contig{$chr2}{rdepth}<0.75 && !exists $p1->{het}) {
#					if (!defined $aa) {
#						$p2->{het} = 1;
#					} else {
#						my $piden = $aa->[M6_PIDEN];
#						if ($piden>97) {$p2->{het} = 2;}
#					}
#				}
#			}
#		}
#	}
#}
#
## -----------------------------------------------
## set cluster size2, NOT count goldfish heterozygous
#my @clust2;
#foreach my $sp (keys(%protein)) {
#	foreach my $pid (keys(%{$protein{$sp}})) {
#		my $aa = $protein{$sp}{$pid};
#		my $tid = $aa->{transcript_id};
#		my $cid = $aa->{clust_id};
#		if (exists $aa->{het}) {next;}
#		push @{$clust2[$cid]{anchor}{$sp}}, [$pid, $tid];
#		$clust2[$cid]{size}++;
#	}
#}
##------------------------------------------------------
#my @clust_ids;
#for (my $cid=0; $cid<@clust2; $cid++) {
#	if (!defined $clust2[$cid]) { next; }
#	my $v = $clust2[$cid]{anchor};
#	my $size = $clust2[$cid]{size};
#	if ($size>2) { push @clust_ids, $cid; }
#}
##------------------------------------------------------

my $sclust_dir="$out.scluster";
if (! -d "$sclust_dir") { mkdir $sclust_dir; }

my @clust_ids;
my @clust2;
for (my $cid=0; $cid<@clust; $cid++) {
	my $c = $clust[$cid];
	if (!defined $c) { next; }
	if ($c->{size}<=2) { next; }
	push @clust_ids, $cid;
}

if (! -d "$clust_dir/0" || ! -d "$sclust_dir/0") {
#-------------------------------------
# Output sequence for each cluster
#-------------------------------------
#{{{
# Load seq
	$in_file = "all5.pep.fasta";
#$in_file = "/data/chenz11/goldfish/11549472/sergey_canu70x/arrow/maker4a/ens85_fish_cc_gc_gf.pep.fasta";
	my $prot_seq = fasta::load_all_seq_hash_ref($in_file, -1); 
	foreach my $sp (keys(%protein)) {
		foreach my $pid (keys(%{$protein{$sp}})) {
			if (exists $prot_seq->{$pid}) {
				$protein{$sp}{$pid}{seq} = $prot_seq->{$pid}{seq};
			}
		}
	}
	$in_file = "all5.rna.fasta";
# $in_file = "/data/chenz11/goldfish/11549472/sergey_canu70x/arrow/maker4a/ens85_fish_cc_gc_gf.rna.fasta";
	my $transcript_seq = fasta::load_all_seq_hash_ref($in_file, -1); 
	foreach my $sp (keys(%transcript)) {
		foreach my $tid (keys(%{$transcript{$sp}})) {
			if (exists $transcript_seq->{$tid}) {
				$transcript{$sp}{$tid}{seq} = $transcript_seq->{$tid}{seq};
			} else {
				delete $transcript{$sp}{$tid};
				my $pid = $transcript{$sp}{$tid}{protein_id};
				if (exists $protein{$sp}{$pid}) {
					delete $protein{$sp}{$pid};
				}
				warn "$tid head no sequence. Delete it.\n";
			}
		}
	}

	if (! -d "$clust_dir/0") {
		for (my $cid=0; $cid<@clust; $cid++) {
			my $c = $clust[$cid];
			if (!defined $c) { next; }
			if ($c->{size}<=2) { next; }
			my $i0 = int($cid / 100);
			my $i1 = $cid % 100;
			my $dir = "$clust_dir/$i0";
			if (! -d $dir) { mkdir $dir; }
			$dir.="/$cid";
			if (! -d $dir) { mkdir $dir; }
			my $out_file1 = "$dir/pep.fa";
			my $out_file2 = "$dir/cdna.fa";
			my $fout1 = base_io::czl_open($out_file1, "w") or die "Fail to create $out_file1\n";
			my $fout2 = base_io::czl_open($out_file2, "w") or die "Fail to create $out_file2\n";
			foreach my $sp (sort keys(%{$c->{anchor}})) {
				foreach my $id (@{$c->{anchor}{$sp}}) {
					my $g = $gene{$sp}{$id};
					my $pid = $gene{$sp}{$id}{protein_id};
					my $tid = $gene{$sp}{$id}{transcript_id};
					if (!defined $pid) {
						print "$tid undef pid\n";
						next;
					} elsif (!defined $tid) {
						print "$pid undef tid\n";
						next;
					}
					my $p = $protein{$sp}{$pid};
					my $t = $transcript{$sp}{$tid};
					if (!defined $p->{seq}) {
						print "$pid undef seq\n";
						next;
					} elsif (!defined $t->{seq}) {
						print "$tid undef seq\n";
						next;
					}
					print $fout1 ">${sp}__$g->{idx}\n$p->{seq}\n";
					print $fout2 ">${sp}__$g->{idx}\n$t->{seq}\n";
				}
			}
			close $fout1;
			close $fout2;
		}
	}

	if (! -d "$sclust_dir/0") {
		for (my $cid=0; $cid<@sclust; $cid++) {
			my $c = $sclust[$cid];
			my $i0 = int($cid / 100);
			my $i1 = $cid % 100;
			my $dir = "$sclust_dir/$i0";
			if (! -d $dir) { mkdir $dir; }
			$dir .= "/$cid";
			if (! -d $dir) { mkdir $dir; }
			my $out_file1 = "$dir/pep.fa";
			my $out_file2 = "$dir/cdna.fa";
			my $fout1 = base_io::czl_open($out_file1, "w") or die "Fail to create $out_file1\n";
			my $fout2 = base_io::czl_open($out_file2, "w") or die "Fail to create $out_file2\n";
			foreach my $sp (sort keys(%{$c->{anchor}})) {
				foreach my $id (@{$c->{anchor}{$sp}}) {
					my $g = $gene{$sp}{$id};
					my $pid = $gene{$sp}{$id}{protein_id};
					my $tid = $gene{$sp}{$id}{transcript_id};
					if (!defined $pid) {
						print "$tid undef pid\n";
						next;
					} elsif (!defined $tid) {
						print "$pid undef tid\n";
						next;
					}
					my $p = $protein{$sp}{$pid};
					my $t = $transcript{$sp}{$tid};
					if (!defined $p->{seq}) {
						print "$pid undef seq\n";
						next;
					} elsif (!defined $t->{seq}) {
						print "$tid undef seq\n";
						next;
					}
					print $fout1 ">${sp}__$g->{idx}\n$p->{seq}\n";
					print $fout2 ">${sp}__$g->{idx}\n$t->{seq}\n";
				}
			}
			close $fout1;
			close $fout2;
		}
	}
# }}}
}


# build trees for each cluster
$out_file="r01.mafft.sh";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
$out_file="r01.clustalo.sh";
$fout1 = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
foreach my $cid (@clust_ids) {
	my $c = $clust[$cid];
	my $i0 = int($cid / 100);
	my $i1 = $cid % 100;
	my $dir = "$clust_dir/$i0/$cid";
#	print $fout "mafft --maxiterate 1000 --localpair $clust_dir/$i0/$cid.pep.fa > $clust_dir/$i0/$cid.pep.mafft.fa;    ";
	print $fout "mafft --retree 1 $dir/pep.fa > $dir/pep.mafft.fa;    ";
#	print $fout "mafft --maxiterate 1000 --localpair $clust_dir/$i0/$cid.cdna.fa > $clust_dir/$i0/$cid.cdna.mafft.fa    ";
	print $fout "mafft --retree 1 $dir/cdna.fa > $dir/cdna.mafft.fa;    ";
	print $fout "tranalign -asequence $dir/cdna.fa -bsequence $dir/pep.mafft.fa -outseq $dir/pep.mafft.codon.fa;    ";
	print $fout "Gblocks $dir/pep.mafft.codon.fa -t=c;    ";
	print $fout "\n";

	print $fout1 "clustalo -i $dir/pep.fa -o $dir/pep.clustalo.fa --full --guidetree-out $dir/pep.clustalo.tree --iterations 10 --force;    ";
	print $fout1 "clustalo -i $dir/cdna.fa -o $dir/cdna.clustalo.fa --full --guidetree-out $dir/cdna.clustalo.tree --iterations 10 --force;    ";
	print $fout1 "tranalign -asequence $dir/cdna.fa -bsequence $dir/pep.clustalo.fa -outseq $dir/pep.clustalo.codon.fa;    ";
	print $fout1 "Gblocks $dir/pep.clustalo.codon.fa -t=c;      ";
	print $fout1 "\n";
}
close $fout;
close $fout1;
# system("swarm --logdir=$log_dir -g 2 -t 1 -m mafft --time=0:30:00 --partition=norm -b 100 -f r01.mafft.sh -J r01.mafft");
# system("swarm --logdir=$log_dir -g 2 -t 1 -m clustalo --time=0:30:00 --partition=norm -b 100 -f r01.clustalo.sh -J r01.clustalo");
$out_file="r01.sclust.mafft.sh";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
$out_file="r01.sclust.clustalo.sh";
$fout1 = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
for (my $cid=0; $cid<@sclust; $cid++) {
	my $c = $sclust[$cid];
	my $i0 = int($cid / 100);
	my $i1 = $cid % 100;
	my $dir = "$sclust_dir/$i0/$cid";
#	print $fout "mafft --maxiterate 1000 --localpair $clust_dir/$i0/$cid.pep.fa > $clust_dir/$i0/$cid.pep.mafft.fa;    ";
	print $fout "mafft --retree 1 $dir/pep.fa > $dir/pep.mafft.fa;    ";
#	print $fout "mafft --maxiterate 1000 --localpair $dir/cdna.fa > $dir/cdna.mafft.fa    ";
	print $fout "mafft --retree 1 $dir/cdna.fa > $dir/cdna.mafft.fa;    ";
	print $fout "tranalign -asequence $dir/cdna.fa -bsequence $dir/pep.mafft.fa -outseq $dir/pep.mafft.codon.fa;    ";
	print $fout "Gblocks $dir/pep.mafft.codon.fa -t=c;   ";
	print $fout "\n";

	print $fout1 "clustalo -i $dir/pep.fa -o $dir/pep.clustalo.fa --full --guidetree-out $dir/pep.clustalo.tree --iterations 10 --force;    ";
	print $fout1 "clustalo -i $dir/cdna.fa -o $dir/cdna.clustalo.fa --full --guidetree-out $dir/cdna.clustalo.tree --iterations 10 --force;    ";
	print $fout1 "tranalign -asequence $dir/cdna.fa -bsequence $dir/pep.clustalo.fa -outseq $dir/pep.clustalo.codon.fa;    ";
	print $fout1 "Gblocks $dir/pep.clustalo.codon.fa -t=c;  ";
	print $fout1 "\n";
}
close $fout;
close $fout1;

# ----------------------------------------
# build trees
# ----------------------------------------
$out_file="r02.mafft.raxml.sh";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
foreach my $cid (@clust_ids) {
	my $c = $clust[$cid];
	my $i0 = int($cid / 100);
	my $i1 = $cid % 100;
	my $dir = "$clust_dir/$i0/$cid";
	print $fout "cd $dir; raxmlHPC -s cdna.mafft.fa -T 1 -p 138927 -n cdna.mafft -m GTRCAT -c 4 ;  cd ../../../;  ";
	print $fout "cd $dir; raxmlHPC -s pep.mafft.fa  -T 1 -p 138927 -n pep.mafft  -m PROTGAMMAAUTO ;  cd ../../../;  ";
	print $fout "\n";
}
close $fout;
$out_file="r02.sclust.mafft.raxml.sh";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
$out_file="r02.sclust.mafft.raxml.og.sh";
$fout1 = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
for (my $cid=0; $cid<@sclust; $cid++) {
	my $c = $sclust[$cid];
	my $ca = $c->{anchor};
	my $i0 = int($cid / 100);
	my $i1 = $cid % 100;
	my $dir = "$sclust_dir/$i0/$cid";
	my $og;
	if (@{$ca->{ENSAMX}} > 0 ) {
		my $g = $gene{ENSAMX}{$ca->{ENSAMX}[0]};
		$og = "ENSAMX__$g->{idx}";
	} else {
		my $g = $gene{ENSDAR}{$ca->{ENSDAR}[0]};
		$og = "ENSDAR__$g->{idx}";
	}
	print $fout "cd $dir; raxmlHPC -s cdna.mafft.fa -T 1 -p 138927 -n cdna.mafft -m GTRCAT -c 4 ; cd ../../../ ;  ";
	print $fout "cd $dir; raxmlHPC -s pep.mafft.fa  -T 1 -p 138927 -n pep.mafft  -m PROTGAMMAAUTO ; cd ../../../ ;  ";
	print $fout "\n";
	print $fout1 "cd $dir; raxmlHPC -s cdna.mafft.fa -T 1 -p 138927 -n cdna.mafft.og -o $og -m GTRCAT -c 4 ; cd ../../../ ;  ";
	print $fout1 "cd $dir; raxmlHPC -s pep.mafft.fa  -T 1 -p 138927 -n pep.mafft.og  -o $og -m PROTGAMMAAUTO ; cd ../../../ ;  ";
	print $fout1 "\n";
}
close $fout;
close $fout1;
# ----------------------------------------

# ----------------------------------------
# compute dn/ds use PAML codeml
# ----------------------------------------
$out_file="r03.sclust.mafft.dnds.sh";
if ($rerun || ! -f $out_file) {
	$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
	for (my $cid=0; $cid<@sclust; $cid++) {
		my $c = $sclust[$cid];
		my $i0 = int($cid / 100);
		my $i1 = $cid % 100;
		my $dir = "$sclust_dir/$i0/$cid";
		my $seq_file = "pep.mafft.codon.pamlphylip";
		my $out_file = "pep.mafft.codon.codemlout";
		my $tree_file= "RAxML_bestTree.pep.mafft";
		# set PAML par: seqfile  outfile  treefile  clock  model  NSite aaDist
		my %par = %paml_par;
		$par{seqfile}  = $seq_file;
		$par{outfile}  = $out_file;
		$par{treefile} = $tree_file;
		$par{runmode}  = -2; # pairwise
		$par{clock}    = 0; # no clock
		$par{model}    = 0; # 1 omega for all branch
		$par{NSsites}  = "0";
		$par{aaDist}   = 0;
		my $out_file1="$dir/codeml.dnds.ctl";
		$fout1 = base_io::czl_open($out_file1, "w") or die "Fail to create $out_file1\n";
		foreach my $u (keys(%par)) {
			if (defined $par{$u}) { print $fout1 "$u = $par{$u}\n"; }
		}
		close $fout1;

		if (! -f "$dir/$seq_file") {
			if (! -f "$dir/pep.mafft.codon.fa") {
				warn "$dir/pep.mafft.codon.fa NOT exists.\n";
				next;
			}
			my $seq = fasta::load_all_seq_ref("$dir/pep.mafft.codon.fa", -1); 
			my $l = length($seq->[0]{seq});
			$out_file1="$dir/$seq_file";
			$fout1 = base_io::czl_open($out_file1, "w") or die "Fail to create $out_file1\n";
			print $fout1 $#{$seq}+1, "\t", $l, "\n";
			foreach my $s (@$seq) {
				$s->{seq} =~ s/(\S\S\S)/$1 /g;
				print $fout1 "$s->{name}\n$s->{seq}\n";
			}
			close $fout1;
		}

		print $fout "cd $dir; codeml codeml.dnds.ctl; cd ../../../  ;   ";
		print $fout "\n";
	}
	close $fout;
}

{
	my %fh;
	my $run = 0;
	my $n=0;
	foreach my $sp1 (@sps) {
		foreach my $sp2 (@sps) {
			if ($sp1 gt $sp2) { next; }
			$out_file="$out.sclust.mafft.dS.$sp1.$sp2.txt";
			$n++;
			if (-f $out_file) { $run++; next; }
			$fh{$sp1}{$sp2} = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
			print {$fh{$sp1}{$sp2}} "#chr1\tID1\tName1\tchr2\tID2\tName2\tN\tS\tdN\tdS\tdN/dS\n";
		}
	}
	if ($run<$n) {
		for (my $cid=0; $cid<@sclust; $cid++) {
			my $c = $sclust[$cid];
			my $ca = $c->{anchor};
			my $i0 = int($cid / 100);
			my $i1 = $cid % 100;
			my $dir = "$sclust_dir/$i0/$cid";
			my $in_file = "$dir/rst";
			if (! -f $in_file) { next; }
			my @ids;
			foreach my $sp (sort keys(%$ca)) {
				foreach my $id (@{$ca->{$sp}}) { push @ids, [$sp, $id]; }
			}
			$fin = base_io::czl_open($in_file, "r") or die "Fail to open $in_file\n";
			my $is_data=0;
			while (<$fin>) {
				if (m/^seq\s+seq/) {
					$is_data=1;
				} elsif ($is_data) {
					s/^\s+//;
					my ($seq1, $seq2, $N, $S, $dN, $dS, $dNdS, @other) = split /\s+/, $_;
					if (!defined $seq1 || $seq1=~m/^\s*$/) { next; }
					if (!defined $seq2 || $seq2=~m/^\s*$/) { next; }
					$seq1--;
					$seq2--;
					my $sp1 = $ids[$seq2][0];
					my $id1 = $ids[$seq2][1];
					my $sp2 = $ids[$seq1][0];
					my $id2 = $ids[$seq1][1];
					if (!exists $fh{$sp1} || !exists $fh{$sp1}{$sp2}) { next; }
					my $g1 = $gene{$sp1}{$id1};
					my $g2 = $gene{$sp2}{$id2};

					print {$fh{$sp1}{$sp2}} "$g1->{seq_name}\t$id1\t$g1->{name}\t$g2->{seq_name}\t$id2\t$g2->{name}\t$N\t$S\t$dN\t$dS\t$dNdS\n";
				}
			}
			close $fin;
		}
	}
	foreach my $sp1 (@sps) {
		foreach my $sp2 (@sps) {
			if ($sp1 gt $sp2) { next; }
			if (exists $fh{$sp1} && exists $fh{$sp1}{$sp2}) { close $fh{$sp1}{$sp2}; }
		}
	}
}


# --------------------------------------------
# Fetch all 1:1:1:2:2 group with same chromosome set
# --------------------------------------------
$n=0;
for (my $scid=0; $scid<@sclust; $scid++) {
	my $sc = $sclust[$scid];
	my $sca= $sc->{anchor};
	if ($sca->{ENSAMX}==0) { next; }
	if ($sca->{CYPCAR}<2) { next; }
	if ($sca->{carAur}<2) { next; }
	my $gf_g0 = $gene{carAur}{$sca->{carAur}[0]};
	my $gf_g1 = $gene{carAur}{$sca->{carAur}[1]};
	my $gf_chr0 = $gf_g0->{seq_name};
	my $gf_chr1 = $gf_g1->{seq_name};
	if ($gf_chr0=~m/LG([0-9]+)/)  { $gf_chr0=$1; } else {next; }
	if ($gf_chr1=~m/LG([0-9]+)/)  { $gf_chr1=$1; } else {next; }
	my $zf_g  = $gene{ENSDAR}{$sca->{ENSDAR}[0]};
	my $zf_chr= $zf_g->{seq_name};
	my $is_match=0;
	for (my $i=0; $i<@{$match_chr{carAur}}; $i++) {
		if ($match_chr{carAur}[$i] eq $gf_chr0) { $is_match++; last; }
	}
	for (my $i=0; $i<@{$match_chr{carAur}}; $i++) {
		if ($match_chr{carAur}[$i] eq $gf_chr1) { $is_match++; last; }
	}
	if ($is_match<2) { next; }

	if ($gf_chr0 > $gf_chr1) {
		@{$sca->{carAur}} = reverse @{$sca->{carAur}};
		($gf_g0, $gf_g1) = ($gf_g1, $gf_g0);
		($gf_chr0, $gf_chr1) = ($gf_chr1, $gf_chr0);
	}
	my $mchr = "chr$zf_chr:$gf_chr0:$gf_chr1";
}
$out_file = "$out.chr_syn.txt";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
close $fout;
exit 0;

# --------------------------------------------
# Fetch all 1:1:1:2:2 group with syn_id set
# --------------------------------------------
my %syn;
my @ctype_count;
foreach my $cid (@clust_ids) {
	my $c = $clust2[$cid];
	my $syn_good=1;
	my @syn_ids;
	my %syn_id_hash;
	my %syn_ids;
	if ($c->{size}>=20) {next;}
	foreach my $sp (qw(ENSDAR CTEIDE CYPCAR carAur)) {
		%{$syn_id_hash{$sp}} = ();
		foreach my $bb (@{$c->{anchor}{$sp}}) {
			my $pid = $bb->[0];
			my $syn_id = $protein{$sp}{$pid}{syn_id};
			$syn_id_hash{$sp}{$syn_id}++;
		}
		@{$syn_ids{$sp}} = sort keys(%{$syn_id_hash{$sp}});
	}
	if (@{$syn_ids{ENSDAR}}>2) { next; }
	if (@{$syn_ids{CTEIDE}}>2) { next; }
	if (@{$syn_ids{CYPCAR}}>6) { next; }
	if (@{$syn_ids{carAur}}>6) { next; }

	my $syn_ids1 = $syn_ids{ENSDAR};
	if (@{$syn_ids1}==0) { @{$syn_ids1} = ("-2"); }

	$syn_ids1 = $syn_ids{CTEIDE};
	if (@{$syn_ids1}==0) { @{$syn_ids1} = ("-2"); }

	$syn_ids1 = $syn_ids{CYPCAR};
	if (@{$syn_ids1}==0) { @{$syn_ids1} = (["-2","-2"]); }
	elsif (@{$syn_ids1}==1) { my $j = $syn_ids1->[0]; @{$syn_ids1} = ([$j,"-2"]); }
	else {
		my @syn_ids2;
		for (my $j1=0; $j1<@$syn_ids1; $j1++) {
			for (my $j2=$j1+1; $j2<@$syn_ids1; $j2++) {
				push @syn_ids2, [ $syn_ids1->[$j1], $syn_ids1->[$j2] ]; 
			}
		}
		@{$syn_ids1} = @syn_ids2;
	}

	$syn_ids1 = $syn_ids{carAur};
	if (@{$syn_ids1}==0) { @{$syn_ids1} = (["-2","-2"]); }
	elsif (@{$syn_ids1}==1) { my $j = $syn_ids1->[0]; @{$syn_ids1} = ([$j,"-2"]); }
	else {
		my @syn_ids2;
		for (my $j1=0; $j1<@$syn_ids1; $j1++) {
			for (my $j2=$j1+1; $j2<@$syn_ids1; $j2++) {
				push @syn_ids2, [ $syn_ids1->[$j1], $syn_ids1->[$j2] ]; 
			}
		}
		@{$syn_ids1} = @syn_ids2;
	}

	for (my $k0=0; $k0<@{$syn_ids{ENSDAR}}; $k0++) {
	for (my $k1=0; $k1<@{$syn_ids{CTEIDE}}; $k1++) {
	for (my $k2=0; $k2<@{$syn_ids{CYPCAR}}; $k2++) {
	for (my $k3=0; $k3<@{$syn_ids{carAur}}; $k3++) {
		my $key = "$syn_ids{ENSDAR}[$k0],$syn_ids{CTEIDE}[$k1]";
		$key .= ",$syn_ids{CYPCAR}[$k2][0],$syn_ids{CYPCAR}[$k2][1]";
		$key .= ",$syn_ids{carAur}[$k3][0],$syn_ids{carAur}[$k3][1]";
		push @{$syn{$key}}, $cid;
	}
	}
	}
	}
}
foreach my $key (keys(%syn)) {
	if (@{$syn{$key}}<=1) { delete $syn{$key}; }
	foreach my $cid (@{$syn{$key}}) {
		push @{$clust2[$cid]{syn_id2}}, $cid;

#		if (exists $clust2[$cid]{syn_id2}) {
#			my $key1 = $clust2[$cid]{syn_id2};
#			my $n1=0;
#			my $n=0;
#			foreach my $syn_id (split ',', $key1) {
#				if ($syn_id>=0) { $n1++; }
#			}
#			foreach my $syn_id (split ',', $key) {
#				if ($syn_id>=0) { $n++; }
#			}
#			if ($n1<$n) {
#				$clust2[$cid]{syn_id2} = $key;
#				$clust[$cid]{syn_id2} = $key;
#			}
#		} else {
#			$clust2[$cid]{syn_id2} = $key;
#			$clust[$cid]{syn_id2} = $key;
#		}
	}
}

#foreach my $cid (@clust_ids) {
#	my $c = $clust2[$cid];
#	if (exists $c->{type} && $c->{type}==CT_SYN1122) { next; }
#	my $type = $c->{type};
#	if (!exists $c->{anchor}{ENSDAR}) { next; }
#	if (!exists $c->{anchor}{ENSAMX}) { next; }
#	if (!exists $c->{anchor}{CTEIDE}) { next; }
#	if (!exists $c->{anchor}{CYPCAR}) { next; }
#	if (!exists $c->{anchor}{carAur}) { next; }
#	my @syn_ids;
#	foreach my $bb0 (@{$c->{anchor}{ENSDAR}}) {
#		my $pid0 = $bb0->[0];
#		my $syn_id0 = $protein{ENSDAR}{$pid0}{syn_id};
#		push @{$syn_ids[0]}, $syn_id0;
#	}
#	foreach my $bb1 (@{$c->{anchor}{CTEIDE}}) {
#		my $pid1 = $bb1->[0];
#		my $syn_id1 = $protein{CTEIDE}{$pid1}{syn_id};
#		push @{$syn_ids[1]}, $syn_id1;
#	}
#	foreach my $bb2 (@{$c->{anchor}{CYPCAR}}) {
#		my $pid20 = $bb2->[0];
#		my $syn_id20 = $protein{CYPCAR}{$pid20}{syn_id};
#		if ($syn_id20 == -1) {next;}
#		foreach my $bb21 (@{$c->{anchor}{CYPCAR}}) {
#			my $pid21 = $bb21->[0];
#			if ($pid20 eq $pid21) { next; }
#			my $syn_id21 = $protein{CYPCAR}{$pid21}{syn_id};
#			if ($syn_id21 == -1) {next;}
#			if ($syn_id20 > $syn_id21) { next; }
#			push @{$syn_ids[2]}, [$syn_id20,$syn_id21];
#		}
#	}
#	foreach my $bb2 (@{$c->{anchor}{carAur}}) {
#		my $pid20 = $bb2->[0];
#		my $syn_id20 = $protein{carAur}{$pid20}{syn_id};
#		if ($syn_id20 == -1) {next;}
#		foreach my $bb21 (@{$c->{anchor}{carAur}}) {
#			my $pid21 = $bb21->[0];
#			if ($pid20 eq $pid21) { next; }
#			my $syn_id21 = $protein{carAur}{$pid21}{syn_id};
#			if ($syn_id21 == -1) {next;}
#			if ($syn_id20 > $syn_id21) { next; }
#			push @{$syn_ids[3]}, [$syn_id20,$syn_id21];
#		}
#	}
#	if (!defined $syn_ids[0] || @{$syn_ids[0]}<1) {next;}
#	if (!defined $syn_ids[1] || @{$syn_ids[1]}<1) {next;}
#	if (!defined $syn_ids[2] || @{$syn_ids[2]}<1) {next;}
#	if (!defined $syn_ids[3] || @{$syn_ids[3]}<1) {next;}
#	for (my $k0=0; $k0<@{$syn_ids[0]}; $k0++) {
#		for (my $k1=0; $k1<@{$syn_ids[1]}; $k1++) {
#			for (my $k2=0; $k2<@{$syn_ids[2]}; $k2++) {
#				for (my $k3=0; $k3<@{$syn_ids[3]}; $k3++) {
#					my $key = $syn_ids[0][$k0];
#					$key .= "|$syn_ids[1][$k1]";
#					$key .= "|$syn_ids[2][$k2][0]|$syn_ids[2][$k2][1]";
#					$key .= "|$syn_ids[3][$k3][0]|$syn_ids[3][$k3][1]";
#					if (exists $syn{$key}) {
#						push @{$syn{$key}}, $cid;
#						$c->{type} = CT_SYN1122a;
#						$ctype_count[CT_SYN1122a]++;
#						last;
#					}
#				}
#				if (exists $c->{type} && $c->{type} == CT_SYN1122a) { last; }
#			}
#			if (exists $c->{type} && $c->{type} == CT_SYN1122a) { last; }
#		}
#		if (exists $c->{type} && $c->{type} == CT_SYN1122a) { last; }
#	}
#}
$out_file = "$out.syn.txt";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
foreach my $key (keys(%syn)) {
	my $m = @{$syn{$key}};
	print $fout "$key\t$m\t", join(",",@{$syn{$key}}), "\n";
}
close $fout;
# --------------------------------------------
# Fetch all 1:1:1:2:2 group with same chromosome set
# --------------------------------------------
my %chr_syn;
my @cchr_type_count;
$n=0;
foreach my $cid (@clust_ids) {
	my $c = $clust2[$cid];
	if ($c->{size}>=20) {next;}
	my %seq_name_hash;
	my %seq_names;
	foreach my $sp (qw(ENSDAR CTEIDE CYPCAR carAur)) {
		%{$seq_name_hash{$sp}} = ();
		foreach my $bb (@{$c->{anchor}{$sp}}) {
			my $pid = $bb->[0];
			my $seq_name = $protein{$sp}{$pid}{seq_name};
			$seq_name_hash{$sp}{$seq_name}++;
		}
		@{$seq_names{$sp}} = sort keys(%{$seq_name_hash{$sp}});
	}
	my $seq_names1 = $seq_names{ENSDAR};
	if (@{$seq_names1}==0) { @{$seq_names1} = ("."); }

	$seq_names1 = $seq_names{CTEIDE};
	if (@{$seq_names1}==0) { @{$seq_names1} = ("."); }

	$seq_names1 = $seq_names{CYPCAR};
	if (@{$seq_names1}==0) { @{$seq_names1} = ([".","."]); }
	elsif (@{$seq_names1}==1) { my $j = $seq_names1->[0]; @{$seq_names1} = ([$j,"."]); }
	else {
		my @seq_names2;
		for (my $j1=0; $j1<@$seq_names1; $j1++) {
			for (my $j2=$j1+1; $j2<@$seq_names1; $j2++) {
				push @seq_names2, [ $seq_names1->[$j1], $seq_names1->[$j2] ]; 
			}
		}
		@{$seq_names1} = @seq_names2;
	}

	$seq_names1 = $seq_names{carAur};
	if (@{$seq_names1}==0) { @{$seq_names1} = ([".","."]); }
	elsif (@{$seq_names1}==1) { my $j = $seq_names1->[0]; @{$seq_names1} = ([$j,"."]); }
	else {
		my @seq_names2;
		for (my $j1=0; $j1<@$seq_names1; $j1++) {
			for (my $j2=$j1+1; $j2<@$seq_names1; $j2++) {
				push @seq_names2, [ $seq_names1->[$j1], $seq_names1->[$j2] ]; 
			}
		}
		@{$seq_names1} = @seq_names2;
	}

	for (my $k0=0; $k0<@{$seq_names{ENSDAR}}; $k0++) {
	for (my $k1=0; $k1<@{$seq_names{CTEIDE}}; $k1++) {
	for (my $k2=0; $k2<@{$seq_names{CYPCAR}}; $k2++) {
	for (my $k3=0; $k3<@{$seq_names{carAur}}; $k3++) {
		my $key = "$seq_names{ENSDAR}[$k0],$seq_names{CTEIDE}[$k1]";
		$key .= ",$seq_names{CYPCAR}[$k2][0],$seq_names{CYPCAR}[$k2][1]";
		$key .= ",$seq_names{carAur}[$k3][0],$seq_names{carAur}[$k3][1]";
		push @{$chr_syn{$key}}, $cid;
	}
	}
	}
	}
}
foreach my $key (keys(%chr_syn)) {
	if (@{$chr_syn{$key}}<=2) {
		delete $chr_syn{$key};
	}
}
$out_file = "$out.chr_syn.txt";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
foreach my $key (keys(%chr_syn)) {
	my $m = @{$chr_syn{$key}};
	print $fout "$key\t$m\t", join(",",@{$chr_syn{$key}}), "\n";
}
close $fout;
# ---------------------------------------------------

my @cluster3;
my @clust_ids3;
{
# ---------------------------------------------------
# filter low identity sequence between goldfish and common carp
# ---------------------------------------------------
# {{{
$in_file="../p.blast/sp5.pairs/CYPCAR.carAur.f3.join.m6.gz";
$fin = base_io::czl_open($in_file, "r") or die "Fail to open $in_file\n";
my %align;
while(<$fin>) {
	chomp;
	my @t = split "\t",$_,-1;
	my ($sp1,$id1) = split /\|/,$t[M6_QID],2;
	my ($sp2,$id2) = split /\|/,$t[M6_TID],2;
	if ($sp1 eq $sp2 && $id1 eq $id2) { next; }
	my $bit   = $t[M6_BIT];
	my $piden = $t[M6_PIDEN];
	my ($cov1,$cov2) = ($t[M6_QCOV],$t[M6_TCOV]);
	if (!exists $protein{$sp1}{$id1}) { next; }
	if (!exists $protein{$sp2}{$id2}) { next; }
	if ($sp1 eq 'CYPCAR' && $sp2 eq 'carAur') {
	} elsif ($sp2 eq 'CYPCAR' && $sp1 eq 'carAur') {
	} else { next; }
	if ($piden<80) { next; }
	my %a = (bit=>$bit,
			score=>$t[M6_SCORE],
			iden=>$piden/100,
			cov1=>$cov1,
			cov2=>$cov2,
			);
	if (!exists $align{$id1}{$id2} || $align{$id1}{$id2}{bit}<$a{bit}) {
		$align{$id1}{$id2} = \%a;
	}
	if (!exists $align{$id2}{$id1} || $align{$id2}{$id1}{bit}<$a{bit}) {
		$a{cov1} = $cov2;
		$a{cov2} = $cov1;
		$align{$id2}{$id1} = \%a;
	}
}
close $fin;
my %best;
foreach my $id1 (keys(%align)) {
	foreach my $id2 (keys(%{$align{$id1}})) {
		my $iden = $align{$id1}{$id2}{iden};
		my $cov1 = $align{$id1}{$id2}{cov1};
		my $cov2 = $align{$id1}{$id2}{cov2};
		my $cov = $cov1>$cov2 ? $cov1 : $cov2;
		if (!exists $best{$id1}) {
			$best{$id1}{id}{$id2}++;
			$best{$id1}{iden_cov} = $iden*$cov;
		} else {
			if ($best{$id1}{iden_cov} < $iden*$cov*0.9) {
				delete $best{$id1}{id};
				$best{$id1}{id}{$id2}++;
			} elsif ($best{$id1}{iden_cov}*0.9 < $iden*$cov) {
				$best{$id1}{id}{$id2}++;
			}
		}
	}
}
foreach my $cid (@clust_ids) {
	my $c = $clust2[$cid];
#	my $syn_id2;
#	my @syn_ids = split ',', $syn_id2, -1;
#	my %syn_id1;
#	my %syn_id2;
#	if ($syn_ids[2]>=0) { $syn_id1{$syn_ids[2]}++; }
#	if ($syn_ids[3]>=0) { $syn_id1{$syn_ids[3]}++; }
#	if ($syn_ids[4]>=0) { $syn_id2{$syn_ids[4]}++; }
#	if ($syn_ids[5]>=0) { $syn_id2{$syn_ids[5]}++; }

	if (!exists $c->{anchor}{ENSDAR} || @{$c->{anchor}{ENSDAR}}!=1) { next; }
	if (!exists $c->{anchor}{CTEIDE} || @{$c->{anchor}{CTEIDE}}!=1) { next; }
	if (!exists $c->{anchor}{ENSAMX} || @{$c->{anchor}{ENSAMX}}==0 || @{$c->{anchor}{ENSAMX}}>2) { next; }
	my @a1;
	my @a2;
	if (!exists $c->{anchor}{CYPCAR}) {@a1=();}
	else { @a1 = @{$c->{anchor}{CYPCAR}}; }
	if (!exists $c->{anchor}{carAur}) {@a2=();}
	else { @a2 = @{$c->{anchor}{carAur}}; }
	my $n1 = @a1;
	my $n2 = @a2;
	if ($n1==0 && $n2==0) { next; }
	my %p1;
	my %p2;
	foreach my $aa (@a1) { $p1{$aa->[0]}=$aa; }
	foreach my $aa (@a2) { $p2{$aa->[0]}=$aa; }
	my $max_iden=0;
	foreach my $aa1 (@a1) {
		my $pid1=$aa1->[0];
		foreach my $aa2 (@a2) {
			my $pid2=$aa2->[0];
			if (exists $align{$pid1}{$pid2}) {
				my $iden = $align{$pid1}{$pid2}{iden};
				if ($max_iden < $iden) { $max_iden = $iden; }
			}
		}
	}
	if ($max_iden<0.8) { next; }
	my $iden_thres = $max_iden-0.1<0.8?0.8:$max_iden-0.1;
	# remove low-identity anchor
	foreach my $aa1 (@a1) {
		my $pid1=$aa1->[0];
		my $max_iden1=0;
		foreach my $aa2 (@a2) {
			my $pid2=$aa2->[0];
			if (exists $align{$pid1}{$pid2}) {
				my $iden = $align{$pid1}{$pid2}{iden};
				if ($max_iden1 < $iden) { $max_iden1 = $iden; }
			}
		}
		if ($max_iden1<$iden_thres) {
#			my $syn_id = $protein{CYPCAR}{$pid1}{syn_id};
#			if (!exists $syn_id1{$syn_id}) {
				delete $p1{$pid1};
				$n1--;
#			}
		}
	}
	@a1 = values(%p1);
	foreach my $aa2 (@a2) {
		my $pid2=$aa2->[0];
		my $max_iden1=0;
		foreach my $aa1 (@a1) {
			my $pid1=$aa1->[0];
			if (exists $align{$pid2}{$pid1}) {
				my $iden = $align{$pid2}{$pid1}{iden};
				if ($max_iden1 < $iden) { $max_iden1 = $iden; }
			}
		}
		if ($max_iden1<$iden_thres) {
			if ($max_iden1<$iden_thres) {
#				my $syn_id = $protein{carAur}{$pid2}{syn_id};
#				if (!exists $syn_id2{$syn_id}) {
					delete $p2{$pid2};
					$n2--;
#				}
			}
		}
	}
	@a2 = values(%p2);

	my $sp1='CYPCAR';
	my $sp2='carAur';
	# add partial by high-identity anchor
	foreach my $aa (@a2) {
		my $pid2 = $aa->[0];
		foreach my $pid1 (keys(%{$align{$pid2}})) {
			if (exists $p1{$pid1}) { next; }
			if ($protein{CYPCAR}{$pid1}{to_dup} ne "-1" && $protein{CYPCAR}{$pid1}{to_dup} ne ".") { next; }
			my $iden = $align{$pid2}{$pid1};
			my $cid1 = $protein{$sp1}{$pid1}{clust_id};
			if ($iden>=$iden_thres && (exists $best{$pid2}{id}{$pid1} && !exists $protein{$sp1}{$pid1}{clust_id3}) ) {
				if (!exists $p1{$pid1}) {
					$p1{$pid1} = [$pid1, $protein{$sp1}{$pid1}{transcript_id}];
					$n1++;
					$protein{$sp1}{$pid1}{clust_id3} = $cid;
				}
			}
		}
	}
	@a1 = values(%p1);
	foreach my $aa (@a1) {
		my $pid1 = $aa->[0];
		foreach my $pid2 (keys(%{$align{$pid1}})) {
			if (exists $p2{$pid2}) { next; }
			if ($protein{carAur}{$pid2}{to_dup} ne "-1" && $protein{carAur}{$pid2}{to_dup} ne ".") { next; }
			my $iden = $align{$pid1}{$pid2};
			if ($iden>=$iden_thres && (exists $best{$pid1}{id}{$pid2} && !exists $protein{$sp2}{$pid2}{clust_id3}) ) {
				if (!exists $p2{$pid2}) {
					$p2{$pid2} = [$pid2, $protein{$sp2}{$pid2}{transcript_id}];
					$n2++;
					$protein{$sp2}{$pid2}{clust_id3} = $cid;
				}
			}
		}
	}
	@a2 = sort values(%p2);

	if (@a1==0 && @a2==0) { next; }

	$cluster3[$cid]{size}=0;
	$cluster3[$cid]{name}=$c->{name};
	my $c3 = $cluster3[$cid];
	foreach my $sp (qw(ENSDAR ENSAMX CTEIDE)) {
		@{$c3->{anchor}{$sp}} = @{$c->{anchor}{$sp}};
	}
	@{$c3->{anchor}{$sp1}} = @a1;
	@{$c3->{anchor}{$sp2}} = @a2;
	foreach my $sp (@sps) {
		foreach my $id (@{$c3->{anchor}{$sp}}) {
			$protein{$sp}{$id->[0]}{clust_id3} = $cid;
			$c3->{size}++;
		}
	}
	push @clust_ids3, $cid;
}
# }}}

# ---------------------------------------------------
# get goldfish expand genes
# ---------------------------------------------------
# {{{
$out_file = "$out.CYPCAR.expand.cluster.txt";
$fout1 = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
$out_file = "$out.carAur.expand.cluster.txt";
$fout2 = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
$out_file = "$out.carAur_CYPCAR.balance.cluster.txt";
my $fout3 = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
foreach my $cid (@clust_ids3) {
	my $c = $cluster3[$cid];
	my $sp1 = 'CYPCAR';
	my $sp2 = 'carAur';
	my $n1 = 0;
	my $n2 = 0;
	if (exists $c->{anchor}{$sp1}) { $n1 = @{$c->{anchor}{$sp1}}; }
	if (exists $c->{anchor}{$sp2}) { $n2 = @{$c->{anchor}{$sp2}}; }
	if ($n1>$n2) {
		print $fout1 "$cid\t$c->{name}\t$n1\t$n2\n";
	} elsif ($n2>$n1) {
		print $fout2 "$cid\t$c->{name}\t$n1\t$n2\n";
	} elsif ($n2==$n1) {
		print $fout3 "$cid\t$c->{name}\t$n1\t$n2\n";
	}
}
close $fout1;
close $fout2;
close $fout3;
# }}}
}
$out_file = "$out.cluster_size3.txt";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
print $fout "#ClusterID\tGeneName\tSize\tSize2\t", join("\t",@sps), "\n";
foreach my $cid (@clust_ids3) {
	my $c = $clust[$cid];
	my $c2= $clust2[$cid];
	my $c3= $cluster3[$cid];
	my $sp = "ENSDAR";
	my $name = $c->{name};
	print $fout "$cid\t$name\t$c->{size}\t$c2->{size}\t$c3->{size}";
	foreach my $sp (@sps) {
		my @m=(0,0,0);
		if (exists $c->{anchor}{$sp} ) { $m[0] = @{ $c->{anchor}{$sp}}; }
		if (exists $c2->{anchor}{$sp}) { $m[1] = @{$c2->{anchor}{$sp}}; }
		if (exists $c3->{anchor}{$sp}) { $m[2] = @{$c3->{anchor}{$sp}}; }
		print $fout "\t", join("\t", @m);
	}
	print $fout "\n";
}
close $fout;

# ---------------------------------------------------
# concat sequnce for each zebrafish chromosome
# ---------------------------------------------------
{
#	my %is_use;
#	foreach my $key (%chr_syn) {
#		my ($zf, $gc, $cc1, $cc2, $gf1, $gf2) = split ".", $key, -1;
#		if ($zf ne ".") { next; }
#		if ($gf1 ne ".") { next; }
#		if ($gf2 ne ".") { next; }
#		foreach my $cid (@{$chr_syn{$key}}) {
#		}
#	}
	my $n = 0;
	my %seq1;
	foreach my $cid (@clust_ids3) {
		my $c3= $cluster3[$cid];
		my $zf_pid = $c3->{anchor}{ENSDAR}[0][0];
		my $gc_pid = $c3->{anchor}{CTEIDE}[0][0];
		my $cf_pid = $c3->{anchor}{ENSAMX}[0][0];
		my $zf_seq_name = $protein{ENSDAR}{$zf_pid}{seq_name};
		my $chr = $zf_seq_name;
		$chr=~ s/^chr//;
		if ($chr !~ /^[0-9]/) { next; }
		my @gf_seq_name = ( "LG$chr", "LG" . ($chr+25) );
		my @cc_seq_name = ( "LG".($chr*2-1), "LG".($chr*2) );
		my $gf_ids = $c3->{anchor}{carAur};
		my $cc_ids = $c3->{anchor}{CYPCAR};
		if (@$gf_ids < 2) {next; }
		if (@$cc_ids < 2) {next; }

		my $gf_pid1;
		my $gf_pid2;
		foreach my $aa (@$gf_ids) {
			my $pid=$aa->[0];
			my $seq_name = $protein{carAur}{$pid}{seq_name};
			if ($seq_name eq $gf_seq_name[0] && !defined $gf_pid1) { $gf_pid1 = $pid; }
			elsif ($seq_name eq $gf_seq_name[1] && !defined $gf_pid2) { $gf_pid2 = $pid; }
		}
		
		my $cc_pid1;
		my $cc_pid2;
		foreach my $aa (@$cc_ids) {
			my $pid=$aa->[0];
			my $seq_name = $protein{CYPCAR}{$pid}{seq_name};
			if ($seq_name eq $cc_seq_name[0] && !defined $cc_pid1) { $cc_pid1 = $pid; }
			elsif ($seq_name eq $cc_seq_name[1] && !defined $cc_pid2) { $cc_pid2 = $pid; }
		}
		
		if (defined $gf_pid1 && defined $gf_pid2 && defined $cc_pid1 && defined $cc_pid2) {
			my $i0 = int($cid/100);
			$in_file = "$clust_dir/$i0/$cid.pep.mafft.codon.fa";
			my $seq2 = fasta::load_all_seq_hash_ref($in_file);
			foreach my $id (keys(%$seq2)) {
				my ($sp,$idx) = split "__", $id, 2;
				if ($sp eq 'ENSDAR' && $protein{$sp}{$zf_pid}{idx}==$idx) { $seq1{ENSDAR}{$zf_seq_name}.=$seq2->{$id}; }
				elsif ($sp eq 'CTEIDE' && $protein{$sp}{$gc_pid}{idx}==$idx) { $seq1{CTEIDE}{$zf_seq_name}.=$seq2->{$id}; }
				elsif ($sp eq 'ENSAMX' && $protein{$sp}{$cf_pid}{idx}==$idx) { $seq1{ENSAMX}{$zf_seq_name}.=$seq2->{$id}; }
				elsif ($sp eq 'CYPCAR' && $protein{$sp}{$cc_pid1}{idx}==$idx) { $seq1{CYPCAR}{$zf_seq_name.'a'}.=$seq2->{$id}; }
				elsif ($sp eq 'CYPCAR' && $protein{$sp}{$cc_pid2}{idx}==$idx) { $seq1{CYPCAR}{$zf_seq_name.'b'}.=$seq2->{$id}; }
				elsif ($sp eq 'carAur' && $protein{$sp}{$gf_pid1}{idx}==$idx) { $seq1{carAur}{$zf_seq_name.'a'}.=$seq2->{$id}; }
				elsif ($sp eq 'carAur' && $protein{$sp}{$gf_pid2}{idx}==$idx) { $seq1{carAur}{$zf_seq_name.'b'}.=$seq2->{$id}; }
			}
			$n++;
		}
	}
	$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
	foreach my $zf_seq_name (keys(%{$seq1{ENSDAR}})) {
		$out_file = "$out.codon.concat.$zf_seq_name.fa";
		foreach my $sp (qw(ENSAMX ENSDAR CTEIDE)) {
			print $fout ">${sp}__$zf_seq_name\n$seq1{$sp}{$zf_seq_name}\n";
		}
		foreach my $sp (qw(CYPCAR carAur)) {
			print $fout ">${sp}__${zf_seq_name}a\n$seq1{$sp}{$zf_seq_name.'a'}\n";
			print $fout ">${sp}__${zf_seq_name}b\n$seq1{$sp}{$zf_seq_name.'b'}\n";
		}
		close $fout;
	}
}


$out_file = "$out.anchor.annot.txt";
$fout = base_io::czl_open($out_file, "w") or die "Fail to create $out_file\n";
foreach my $sp (keys(%protein)) {
	foreach my $pid (keys(%{$protein{$sp}})) {
		my $aa = $protein{$sp}{$pid};
		print $fout "$sp\t$aa->{seq_name}\t$aa->{begin}\t$aa->{end}\t$aa->{strand}\t$aa->{idx}\t$pid\t$aa->{transcript_id}\t$aa->{gene_id}\t$aa->{to_dup}\t$aa->{pos_idx}\t$aa->{clust_id0}\t$aa->{clust_id}\t$aa->{syn_id}\t";
		if (!exists $aa->{syn_id2}) { print $fout "."; }
		else { print $fout join(";",$aa->{syn_id2}); }
		print $fout "\n";
	}
}
close $fout;
