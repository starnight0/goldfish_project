#####################################
# 1. nearby (distance < 10 exons) GF genes is group into split-gene if they map 
#    to a same ZF gene. Compute the total covered GF and ZF exon lenght and 
#    the fraction each split-gene covered
# 2. GF genes map to the same ZF gene is group in to a tuple
#####################################
use strict;
use warnings;
use Data::Dumper;

my $out = $ARGV[0];
my %tuple;
my %to_zf;
my %bad;
my $sp = "GF";
my %pairs_exon;
my $n=0;
my $m=0;
my $zf_gid0;
my %cur_gid1;
my %cur_exon1;
my %cur_exon2;
my %cur_len1;
my %cur_len2;
my @v;
while(<STDIN>) {
	my @t = split /\t/;
	if ($t[17] eq ".") { next; }
	push @v, \@t; 
}
$n = @v;

my @exon;
my @gene;
for (my $i1=0; $i1<@v; $i1++) {
	my @id = ($v[$i1][3],$v[$i1][17]);
	for (my $k=0; $k<2; $k++) {
		my @id1=split /:/, $id[$k];
		my $gid1=$id1[4];
		my $eid1 = $id[$k];
		if (!exists $exon[$k]{$eid1}) {
			$exon[$k]{$eid1}{gene_id} = $gid1;
			my $l  = $id1[2]-$id1[1];
			$exon[$k]{$eid1}{size} = $l;
			$exon[$k]{$eid1}{map_num}++;
			$gene[$k]{$gid1}{size} += $l;
		}
	}
}

# set uniq mapped position
my @is_uniq = (0) x $n;
my @u;
my $i0=0;
my $i1=0;
for ($i1=1; $i1<@v; $i1++) {
	if ($v[$i1][0] ne $v[$i0][0] || $v[$i1][1] != $v[$i0][1]) {
		if ($i1-$i0==1) { $is_uniq[$i0]=1; }
		else {push @u, [$i0,$i1];}
		$i0=$i1;
	}
}
if ($i1-$i0==1) { $is_uniq[$i0]=1; }
else {push @u, [$i0,$i1];}
##

# set map genes for unique mapped pos0
for (my $i0=0; $i0<@v; $i0++) {
	if (!$is_uniq[$i0]) { next; }
	my @t = @{$v[$i0]};
	my $id1s=$t[3];
   	my $id2s=$t[17];
	my @id1=split /:/, $id1s;
	my @id2=split /:/, $id2s;
	my $gid1=$id1[4];
	my $gid2=$id2[4];
	if (!exists $gene[0]{$gid1}{map_gene}{$gid2}{exons2}{$id2s}) {
		$gene[0]{$gid1}{map_gene}{$gid2}{size2}+=$exon[1]{$id2s}{size};
	}
	$gene[0]{$gid1}{map_gene}{$gid2}{exons2}{$id2s}++;
}

# choose best mapped pos1 if not unique
for (my $ui=0; $ui<@u; $ui++) {
	my $i0 = $u[$ui][0];
	my $i1 = $u[$ui][1];
	my %u;
	my ($gid1,$gid2);
	for (my $i=$i0; $i<$i1; $i++) {
		my @t = @{$v[$i]};
		my $id1s=$t[3];
		my $id2s=$t[17];
		my @id1=split /:/, $id1s;
		my @id2=split /:/, $id2s;
		$gid1=$id1[4];
		$gid2=$id2[4];
		if (exists $gene[0]{$gid1}) {
			if (exists $gene[0]{$gid1}{map_gene}{$gid2}) {
				my $l = $gene[0]{$gid1}{map_gene}{$gid2}{size2};
				$u{$gid2} = [$i, $l];
			}
		}
	}
	my @u1;
	map { push @u1, [$_,@{$u{$_}}]} keys(%u);
	if (@u1==0) {
		for (my $i=$i0; $i<$i1; $i++) { $is_uniq[$i] = 5; }
	} elsif (@u1==1) {
		my $i = $u1[0][1];
		$is_uniq[$i] = 2;
		my $gid2 = $u1[0][0];
		$gene[0]{$gid1}{map_gene}{$gid2}{size2} += $u1[0][2];
	} else {
		@u1 = sort {-($a->[2]<=>$b->[2])} @u1;
		my $max = $u1[0][2];
		if ($u1[1][2]<$max*0.8) {
			my $i = $u1[0][1];
			$is_uniq[$i] = 3;
			my $gid2 = $u1[0][0];
			$gene[0]{$gid1}{map_gene}{$gid2}{size2} += $u1[0][2];
		} else {
			for (my $j=0; $j<@u1; $j++) {
				if ($u1[$j][2]<$max*0.8) { last; }
				my $i = $u1[$j][1];
				$is_uniq[$i] = 4;
				my $gid2 = $u1[0][0];
				$gene[0]{$gid1}{map_gene}{$gid2}{size2} += $u1[0][2];
			}
		}
	}
}
##

my $dist=10;
for (my $i0=0; $i0<@v; $i0++) {
	if (!$is_uniq[$i0]) { next; }
	my @t = @{$v[$i0]};
	my $id1=$t[3];    my $id2=$t[17];
	my @id1=split /:/, $id1;  my @id2=split /:/, $id2;
	my $gid1=$id1[4];  my $gid2=$id2[4];
	$zf_gid0 = $gid2;
	if ($id2 eq ".") { next; }
	my $d = 0;
	my $i2;
	my $i1;
	for ($i1=$i0+1; $i1<@v; $i1++) {
		if (!$is_uniq[$i1]) { next; }
		my $t1 = $v[$i1];
		my $id3=$t1->[3];    my $id4=$t1->[17];
		my @id3=split /:/, $id3;  my @id4=split /:/, $id4;
		my $gid3=$id3[4];  my $gid4=$id4[4];
		if ($gid4 eq $zf_gid0) { $d=0; next; }
		else {
			if (!defined $i2) { $i2=$i1; }
			if ($d>$dist) { $i1-=$d; last; }
			else { $d++; }
		}
	}
	$id1[1]--;   $id2[1]--;
	%cur_len1=();   %cur_gid1=();   %cur_exon1=();
	%cur_len2=();   %cur_exon2=(); 
	for (my $i=$i0; $i<$i1; $i++) {
		if (!$is_uniq[$i]) { next; }
		my $t1 = $v[$i];
		my $id1=$t1->[3];    my $id2=$t1->[17];
		my @id1=split /:/, $id1;  my @id2=split /:/, $id2;
		my $gid1=$id1[4];  my $gid2=$id2[4];
		if ($gid2 ne $zf_gid0) { next; }

		$cur_gid1{$gid1}[2]++; # counts
		$cur_gid1{$gid1}[3]+=$t1->[8]-$t1->[7]; # length
		if (!defined $cur_gid1{$gid1}[1]) {
			$cur_gid1{$gid1}[0] = $id1[0]; # chr
			$cur_gid1{$gid1}[1] = $id1[1]; # begin
		} else {
			if ($cur_gid1{$gid1}[1] > $id1[1]) { $cur_gid1{$gid1}[1] = $id1[1]; }
		}

		if (!exists $cur_exon1{$id1}) { $cur_exon1{$id1}++;  $cur_len1{$gid1}+=$id1[2]-$id1[1]; }
		if (!exists $cur_exon2{$id2}) { $cur_exon2{$id2}++;  $cur_len2{$gid2}+=$id2[2]-$id2[1]; }
	}
	my $l=0;   foreach my $gid1 (keys(%cur_gid1)) { $l+=$cur_len1{$gid1}; }
	my $l2=0;  foreach my $gid2 (keys(%cur_len2)) { $l2+=$cur_len2{$gid2}; }
	my %gs;
	foreach my $gid1 (sort {$cur_gid1{$a}[1]<=>$cur_gid1{$b}[1]} keys(%cur_gid1)) {
		my $a = $cur_gid1{$gid1};
		if ($a->[3]<0.1*$l) { next; }
		$gs{$gid1} = [@$a, sprintf("%.2f", $a->[3]/$l) ];
	}
	if (%gs) { push @{$tuple{$zf_gid0}{$sp}}, [\%gs, $l, $l2]; }

	if (defined $i2) { $i0 = $i2-1; }
}

# remove short tuples
foreach my $gid2 (sort keys(%tuple)) {
	if (!exists $tuple{$gid2}{GF}) { next; }
	my $m = @{$tuple{$gid2}{GF}};
	if ($m<=1) { next; }
	my $l2 = 0;
	my @t = @{$tuple{$gid2}{GF}};
	foreach my $aa (@t) {
		if ($aa->[2]>$l2) { $l2=$aa->[2]; } # total lengh by ZF exons
	}
	my $i0=0;
	for (my $i=0; $i<@t; $i++) {
		my $aa = $t[$i];
		if ($aa->[2]<0.1*$l2) { next; }
		if ($i!=$i0) { $t[$i0] = $t[$i]; }
		$i0++;
	}
	$#{$tuple{$gid2}{GF}} = $i0-1;

}

$n = keys(%tuple);
print STDERR "$sp: $n\n";

# output all ZF,GF tuple
open OUT1, ">${out}tuple";
open OUT2, ">${out}tuple.gene_ids";
foreach my $gid2 (sort keys(%tuple)) {
	if (!exists $tuple{$gid2}{GF}) { next; }
	my $m = @{$tuple{$gid2}{GF}};
	print OUT1 $gid2;
	print OUT2 $gid2;
	for (my $i=0; $i<$m; $i++) {
		print OUT1 "\t";
		print OUT2 "\t";
		my $a1 = $tuple{$gid2}{GF}[$i][0];
		my @u;
		foreach my $gid1 (keys(%$a1)) {
			push @u, join(":",(@{$a1->{$gid1}},$gid1));
		}
		print OUT1 join(',', @u);
		print OUT2 join(',', keys(%$a1));
	}
	print OUT1 "\n";
	print OUT2 "\n";
}
close OUT1;
close OUT2;

open OUT1, ">${out}sp1_sp2_gene_pairs";
open OUT2, ">${out}sp1_sp2_gene_pairs2";
foreach my $gid2 (sort keys(%tuple)) {
	if (!exists $tuple{$gid2}{GF}) { next; }
	my $m = @{$tuple{$gid2}{GF}};
	for (my $i=0; $i<$m; $i++) {
		my $a1 = $tuple{$gid2}{GF}[$i][0];
		print OUT1 $gid2,"\t",join(',', keys(%$a1)), "\n";
		foreach my $gid1 (keys(%$a1)) {
			print OUT2 $gid2,"\t",$gid1,"\n";
		}
	}
}
close OUT1;
close OUT2;

open OUT1, ">${out}sp2_sp2_gene_pairs";
open OUT2, ">${out}sp2_sp2_gene_pairs2";
foreach my $gid2 (sort keys(%tuple)) {
	if (!exists $tuple{$gid2}{GF}) { next; }
	my $m = @{$tuple{$gid2}{GF}};
	if ($m<2) { next; }
	for (my $i=0; $i<$m-1; $i++) {
		my $a1 = $tuple{$gid2}{GF}[$i][0];
		for (my $j=$i+1; $j<$m; $j++) {
			my $a2 = $tuple{$gid2}{GF}[$j][0];
			print OUT1 $gid2,"\t",join(',', keys(%$a1)),"\t",join(',', keys(%$a2)), "\n";
			foreach my $gid1 (keys(%$a1)) {
				foreach my $gid3 (keys(%$a2)) {
					print OUT2 $gid2,"\t",$gid1,"\t",$gid3,"\n";
				}
			}
		}
	}
}
close OUT1;
close OUT2;


#if ( -f "${out}1.txt") { die "OUTFILE exists.\n" }
#if ( -f "${out}2.txt") { die "OUTFILE exists.\n" }
#if ( -f "${out}2.txt") { die "OUTFILE exists.\n" }
my @fh;
my @fh2;
open $fh[1], ">${out}1.txt" or die;
open $fh[2], ">${out}2.txt" or die;
open $fh[3], ">${out}3.txt" or die;
open $fh2[1], ">${out}1.2.txt" or die;
open $fh2[2], ">${out}2.2.txt" or die;
open $fh2[3], ">${out}3.2.txt" or die;
open OUT, ">${out}2.3.txt" or die;
my %pass;
my %tuple_hist;
foreach my $gid2 (sort keys(%tuple)) {
	if (!exists $tuple{$gid2}{GF}) { next; }
	my $m = @{$tuple{$gid2}{GF}};
	$tuple_hist{$m}++;
	my $t = $tuple{$gid2}{GF};
	if ($m>=1 && $m<=3) {
		my $bad=0;
		for (my $i=0; $i<$m; $i++) {
			my $a1 = $tuple{$gid2}{GF}[$i][0];
			foreach my $gid1 (keys(%$a1)) { if ($pass{$gid1}) {$bad=1; last; } }
		}
		if ($bad) { next; }
		print {$fh[$m]} $gid2;
		print {$fh2[$m]} $gid2;
		my $is_out3=1;
		if ($m!=2) { $is_out3=0; }
		else {
			my $max_l1;
			my $max_l2;
			for (my $i=0; $i<$m; $i++) {
				my $l1 = $t->[$i][1];
				my $l2 = $t->[$i][2];
				if (!defined $max_l1 || $max_l1<$l1) { $max_l1 = $l1; }
				if (!defined $max_l2 || $max_l2<$l2) { $max_l2 = $l2; }
			}
			@$t = sort { -($a->[1] <=> $b->[1]) } @$t;
			for (my $i=0; $i<$m; $i++) {
				my $a1 = $t->[$i][0];
				if (keys(%$a1)>1) { $is_out3=0; last; }
#				my @gid1s = sort {-($a1->{$a}[2]<=>$a1->{$b}[2])} keys(%$a1);
#				my $gid10 = $gid1s[0];
#				if ($a1->{$gid10}[3]<$max_l1*0.5) { $is_out3=0; last; }
			}
		}
		if ($is_out3) { print OUT  $gid2; }
		for (my $i=0; $i<$m; $i++) {
			print {$fh[$m]} "\t";
			print {$fh2[$m]} "\t";
			if ($is_out3) { print OUT  "\t"; }
			my $a1 = $tuple{$gid2}{GF}[$i][0];
			my @gid1s = sort {-($a1->{$a}[1]<=>$a1->{$b}[1])} keys(%$a1);
			my $gid10 = $gid1s[0];
			my $j=0;
			foreach my $gid1 (sort {$a1->{$a}[1]<=>$a1->{$b}[1]} keys(%$a1)) {
				if ($j>0) { print {$fh[$m]} ","; }
				$j++;
				print {$fh[$m]} join(":",(@{$a1->{$gid1}},$gid1)); $pass{$gid1}++;
			}
			print {$fh2[$m]} join(":",(@{$a1->{$gid10}},$gid10));
			if ($is_out3) { print OUT join(":",(@{$a1->{$gid10}},$gid10)); }
		} 
		print {$fh[$m]} "\n"; print {$fh2[$m]} "\n";
		if ($is_out3) { print OUT "\n"; }
	}
}
close OUT;
for (my $i=1; $i<=3; $i++) { close $fh[$i]; close $fh2[$i]; }

print STDERR "histogram:\n";
foreach my $m (sort {$a<=>$b} keys(%tuple_hist)) {
	print STDERR "$m\t$tuple_hist{$m}\n";
}
