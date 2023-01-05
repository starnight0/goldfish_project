wd=big/WGD2
mkdir $wd/input
for sp in ZF CC GF; do mkdir $wd/input/$sp; done
mkdir $wd/output


chainFilter -qMinSize=20000 -tMinSize=20000 big/WGD/GF.GF.net.chain | perl -ne '
if (m/^#/) { next; } chomp; 
if (m/^chain/) {
	    my @t=split / /;
		    if ($t[2] eq $t[7]) { $filt=1;}
			    else {$filt=0; print $_, "\n";}
} elsif (!$filt) {
	    print $_, "\n";
}' > big/WGD/GF.GF.net.long20000.chain

chainFilter -qMinSize=20000 -tMinSize=20000 big/WGD/GF.ZF.net.chain > big/WGD2/input/GF.ZF.net.long20000.chain 
#perl -I ~/my_program3/src/perl_pm ~/my_program3/src/utility/czl_chain_remove_TQovl.pl -m T -ovlf 0.2 -i big/WGD2/input/GF.ZF.net.long.chain -o big/WGD2/input/GF.ZF.net.long.non_Tovl.chain
#cat big/WGD2/input/GF.ZF.net.long.non_Tovl.chain | grep chain | awk '{print $3"\t"$6"\t"$7}' | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk '{l+=$3-$2;} END{print l}'
#chainSort big/WGD2/input/GF.ZF.net.long.non_Tovl.chain big/WGD2/input/GF.ZF.net.long.non_Tovl.sorted_by_score.chain 

chainFilter -qMinSize=20000 -tMinSize=20000 big/WGD/CC.ZF.net.chain > big/WGD2/input/CC.ZF.net.long20000.chain 
#perl -I ~/my_program3/src/perl_pm ~/my_program3/src/utility/czl_chain_remove_TQovl.pl -m T -ovlf 0.2 -i big/WGD2/input/CC.ZF.net.long.chain -o big/WGD2/input/CC.ZF.net.long.non_Tovl.chain
#chainSort big/WGD2/input/CC.ZF.net.long.non_Tovl.chain big/WGD2/input/CC.ZF.net.long.non_Tovl.sorted_by_score.chain 

# get exons region for each species
# {{{
cat /home/chenz11/data/zebrafish/ensembl85/Danio_rerio.GRCz10.85.noM.longest_transcript.gtf | perl -ne '
if (m/^#/ || m/^\s*$/) { next; }
s/\s+$//;
my @t = split /\t/;
my %info;
foreach my $a (split /;/, $t[8]) {
	$a=~s/^\s+//; $a=~s/\s+$//;
	my ($u,$v) = split / /, $a;
	$u=~s/^\s+//; $u=~s/\s+$//;
	$v=~s/^\s+//; $v=~s/\s+$//;
	$v=~s/^"(.*)"$/\1/;
	$info{$u} = $v;
}
if ($t[2] eq "gene") {
	$gene_id = $info{gene_id};
} elsif ($t[2] eq "transcript") {
	$transcript_id = $info{transcript_id};
} elsif ($t[2] eq "exon") {
	$exon_id = $info{exon_id};
	$chr = $t[0];
	if ($chr=~m/^[0-9]/) { $chr="chr$chr"; }
	elsif ($chr=~m/^MT$/) { $chr="chrM"; }
	elsif ($chr=~m/^KN/) { $chr="chrUn_$chr"; $chr=~s/\./v/; }
	print join("\t", ($chr, $t[3]-1, $t[4], "$chr:$t[3]:$t[4]:$t[6]:$gene_id:$transcript_id:$exon_id", 0, $t[6])), "\n";
}
' | sort -k1,1 -k2,2n > $wd/input/ZF/exon.bed;

cat big/carAur03.noM.gene.unmasked.gff | perl -ne '
if (m/^#/ || m/^\s*$/) { next; }   s/\s+$//;
my @t = split /\t/;      my %info;
foreach my $a (split /;/, $t[8]) {
	$a=~s/^\s+//; $a=~s/\s+$//;
	my ($u,$v) = split /=/, $a;
	$u=~s/^\s+//; $u=~s/\s+$//;
	$v=~s/^\s+//; $v=~s/\s+$//;
	$v=~s/^"(.*)"$/\1/;
	$info{$u} = $v;
}
if ($t[2] eq "gene") {
	$gene_id = $info{ID};
} elsif ($t[2] eq "transcript" || $t[2]=~m/RNA/) {
	$transcript_id = $info{ID};
	$transcript_id =~ s/\.[0-9]+$//;
} elsif ($t[2] eq "exon") {
	$exon_id = $info{ID};
	$exon_id =~ s/:/__/g;
	$chr = $t[0];
	print join("\t", ($chr, $t[3]-1, $t[4], "$chr:$t[3]:$t[4]:$t[6]:$gene_id:$transcript_id:$exon_id", 0, $t[6])), "\n";
}
' | sort -k1,1 -k2,2n > $wd/input/GF/exon.bed;

cat ~/data/common_carp/NCBI_Cyprinus_carpio/ref_common_carp_genome_top_level.rename_chr.re_id1.noM.gff3 | perl -ne '
if (m/^#/ || m/^\s*$/) { next; }   s/\s+$//;
my @t = split /\t/;      my %info;
foreach my $a (split /;/, $t[8]) {
	$a=~s/^\s+//; $a=~s/\s+$//;
	my ($u,$v) = split /=/, $a;
	$u=~s/^\s+//; $u=~s/\s+$//;
	$v=~s/^\s+//; $v=~s/\s+$//;
	$v=~s/^"(.*)"$/\1/;
	$info{$u} = $v;
}
if ($t[2] eq "gene") {
	$gene_id = $info{ID};
} elsif ($t[2] eq "transcript") {
	$transcript_id = $info{ID};
	$transcript_id =~ s/\.[0-9]+$//;
} elsif ($t[2] eq "exon") {
	$exon_id = $info{ID};
	$chr = $t[0];
	print join("\t", ($chr, $t[3]-1, $t[4], "$chr:$t[3]:$t[4]:$t[6]:$gene_id:$transcript_id:$exon_id", 0, $t[6])), "\n";
}
' | sort -k1,1 -k2,2n > $wd/input/CC/exon.bed;
# }}}

# get CNE region for each species
# {{{
cat big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.GF.notR.bed | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp1;
bedtools subtract -nonamecheck -A -a big/phastcons/run2.net_roast_multic/ZF_GF.GF.notRE.bed -b big/WGD2/tmp/tmp1 | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp2;
cat big/WGD2/tmp/tmp1 big/WGD2/tmp/tmp2 | bedSort stdin big/WGD2/tmp/tmp3; mv big/WGD2/tmp/tmp3  big/WGD2/tmp/tmp1;
bedtools subtract -nonamecheck -A -a big/WGD2/tmp/tmp1 -b big/WGD2/input/GF/exon.bed | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp2;
cat big/WGD2/tmp/tmp2 | awk -F$'\t' -v OFS=$'\t' '{$4="CNE1:"$1":"$2+1":"$3":"$4; $5=0; print}' | grep -v 'chrM' | bedSort stdin big/WGD2/input/GF/CNE.bed

cat big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.CC.notR.bed | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp1;
bedtools subtract -nonamecheck -A -a big/phastcons/run2.net_roast_multic/CC_GF.CC.notR.bed -b big/WGD2/tmp/tmp1 | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp2;
cat big/WGD2/tmp/tmp1 big/WGD2/tmp/tmp2 | bedSort stdin big/WGD2/tmp/tmp3; mv big/WGD2/tmp/tmp3  big/WGD2/tmp/tmp1;
bedtools subtract -nonamecheck -A -a big/WGD2/tmp/tmp1 -b big/WGD2/input/CC/exon.bed | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp2;
cat big/WGD2/tmp/tmp2 | awk -F$'\t' -v OFS=$'\t' '{$4="CNE1:"$1":"$2+1":"$3":"$4; $5=0; print}' | grep -v 'chrM' | bedSort stdin big/WGD2/input/CC/CNE.bed

cat big/phastcons/run2.net_roast_multic/ZF_GC_CC_GF.ZF.notR.bed | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp1;
bedtools subtract -nonamecheck -A -a big/phastcons/run2.net_roast_multic/ZF_GF.ZF.notR.bed -b big/WGD2/tmp/tmp1 | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp2;
cat big/WGD2/tmp/tmp1 big/WGD2/tmp/tmp2 | bedSort stdin big/WGD2/tmp/tmp3; mv big/WGD2/tmp/tmp3  big/WGD2/tmp/tmp1;
bedtools subtract -nonamecheck -A -a big/WGD2/tmp/tmp1 -b big/WGD2/input/ZF/exon.bed | sort -k1,1 -k2,2n > big/WGD2/tmp/tmp2;
cat big/WGD2/tmp/tmp2 | awk -F$'\t' -v OFS=$'\t' '{$4="CNE1:"$1":"$2+1":"$3":"$4; $5=0; print}' | grep -v 'chrM' | bedSort stdin big/WGD2/input/ZF/CNE.bed
# }}}

for sp in CC GF ZF
do
cat big/WGD2/input/$sp/exon.bed | awk -F$'\t' -v OFS=$'\t' '{print $1,$2,$3,"exon:"$4,$5,$6}' > a;
cat a big/WGD2/input/$sp/CNE.bed | sort -k1,1 -k2,2n > big/WGD2/input/$sp/exon_CNE.bed
done

for sp in CC GF
do
crossmap bed big/WGD2/input/$sp.ZF.net.long20000.chain $wd/input/$sp/exon.bed $wd/input/$sp/exon.map.bed; 
cat $wd/input/$sp/exon.map.bed | perl ~/my_program3/src/utility/czl_bed_merge_only_same_id.pl 201 > $wd/input/$sp/exon.map1.bed ;
#cat $wd/input/$sp/exon.map1.bed | perl -ne 'chomp; my @t=split /\t/; my @id=split /:/, $t[3];
#$id[1]--;
#my $l0=$id[2]-$id[1];
#my $l1=$t[2]-$t[1];
#if ($l1>=0.5*$l0) {print $_, "\n"; }
#' > $wd/input/$sp/exon.map1_C50.bed ;
#cat $wd/input/$sp/exon.map1.bed | sort -k1,1 -k2,2n > $wd/input/$sp/exon.map1.sorted.bed;
bedtools intersect -nonamecheck -wao -f 0.5 -a $wd/input/$sp/exon.map1.bed -b $wd/input/ZF/exon.bed > $wd/input/$sp/exon.map1_ZF.bed;
cat $wd/input/$sp/exon.map1_ZF.bed | perl -ne '
chomp; @t=split /\t/;  my $id=$t[3]; $t[3]=".";  my @id=split /:/, $id;
$id[1]--;
print join("\t",(@id[0..2],$id,0,$id[3], @t)), "\n";
' | sort -k1,1 -k2,2n -k13,13 -k14,14n  > $wd/input/$sp/exon.map1_ZF3.bed;
done

# fetch double-retained genes for both CC and GF (quintuple)
perl -e '
my %tuple;
my %to_zf;
my %bad;
foreach my $sp (("CC","GF")) {
	open IN, "<'$wd'/input/$sp/exon.map_C50_ZF.bed" or die;
	my $n=0; my $m=0;
	while(<IN>) {
		my @t = split /\t/;
		if ($t[9] eq ".") { next; }
		my @id1=split /:/, $t[3];  my @id2=split /:/, $t[9];
		my $gid1=$id1[4];  my $gid2=$id2[4];
		if (exists $bad{$gid2}) { next; }
		if (exists $to_zf{$gid1}) {
			if ($to_zf{$gid1} eq $gid2) { next; }
			else {$bad{$gid1}++; $bad{$gid2}++; delete $tuple{$to_zf{$gid1}}; delete $tuple{$gid2}; next;}
		}
		$to_zf{$gid1} = $gid2;

#cat big/WGD/CNE.CC_GF.bed big/WGD/CNE.ZF_GF.bed | bedtools subtract -nonamecheck -A -a big/WGD/CNE.GF_GF.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.GF_GF_only.bed
			cat big/WGD/CNE.ZF_GF.bed big/WGD/CNE.GC_GF.bed | bedtools subtract -nonamecheck -A -a big/WGD/CNE.CC_GF.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.CC_GF.1.bed
			cat big/WGD/CNE.CC_GF.1.bed big/WGD/CNE.ZF_GF.bed | bedtools subtract -nonamecheck -A -a big/WGD/CNE.GC_GF.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.GC_GF.1.bed
			cat big/WGD/CNE.CC_GF.1.bed big/WGD/CNE.GC_GF.1.bed | bedtools subtract -nonamecheck -A -a big/WGD/CNE.ZF_GF.bed -b stdin | sort -k1,1 -k2,2n > big/WGD/CNE.ZF_GF.1.bed
# }}}

for sp in CC GF
do
crossmap bed big/WGD2/input/$sp.ZF.net.long20000.chain $wd/input/$sp/exon.bed $wd/input/$sp/exon.map.bed; 
cat $wd/input/$sp/exon.map.bed | perl ~/my_program3/src/utility/czl_bed_merge_only_same_id.pl 201 > $wd/input/$sp/exon.map1.bed ;
#cat $wd/input/$sp/exon.map1.bed | perl -ne 'chomp; my @t=split /\t/; my @id=split /:/, $t[3];
#$id[1]--;
#my $l0=$id[2]-$id[1];
#my $l1=$t[2]-$t[1];
#if ($l1>=0.5*$l0) {print $_, "\n"; }
#' > $wd/input/$sp/exon.map1_C50.bed ;
#cat $wd/input/$sp/exon.map1.bed | sort -k1,1 -k2,2n > $wd/input/$sp/exon.map1.sorted.bed;
bedtools intersect -nonamecheck -wao -f 0.5 -a $wd/input/$sp/exon.map1.bed -b $wd/input/ZF/exon.bed > $wd/input/$sp/exon.map1_ZF.bed;
cat $wd/input/$sp/exon.map1_ZF.bed | perl -ne '
chomp; @t=split /\t/;  my $id=$t[3]; $t[3]=".";  my @id=split /:/, $id;
$id[1]--;
print join("\t",(@id[0..2],$id,0,$id[3], @t)), "\n";
' | sort -k1,1 -k2,2n -k13,13 -k14,14n  > $wd/input/$sp/exon.map1_ZF3.bed;
done

# fetch double-retained genes for both CC and GF (quintuple)
perl -e '
my %tuple;
my %to_zf;
my %bad;
foreach my $sp (("CC","GF")) {
	open IN, "<'$wd'/input/$sp/exon.map_C50_ZF.bed" or die;
	my $n=0; my $m=0;
	while(<IN>) {
		my @t = split /\t/;
		if ($t[9] eq ".") { next; }
		my @id1=split /:/, $t[3];  my @id2=split /:/, $t[9];
		my $gid1=$id1[4];  my $gid2=$id2[4];
		if (exists $bad{$gid2}) { next; }
		if (exists $to_zf{$gid1}) {
			if ($to_zf{$gid1} eq $gid2) { next; }
			else {$bad{$gid1}++; $bad{$gid2}++; delete $tuple{$to_zf{$gid1}}; delete $tuple{$gid2}; next;}
		}
		$to_zf{$gid1} = $gid2;
		push @{$tuple{$gid2}{$sp}}, $gid1;
	}
	close IN;
	my $n = keys(%tuple);
	print STDERR "$sp: $n\n";
}
foreach my $gid2 (sort keys(%tuple)) {
	if (exists $tuple{$gid2}{CC} && @{$tuple{$gid2}{CC}}==2 
			&& exists $tuple{$gid2}{GF} && @{$tuple{$gid2}{GF}}==2) {
		print join("\t", ($gid2,@{$tuple{$gid2}{CC}}, @{$tuple{$gid2}{GF}})), "\n";
	}
}
' > $wd/input/quintuple.txt


# fetch double-retained genes for GF (triplets)
for sp in CC GF
do
	cat $wd/input/$sp/exon.map1_ZF3.bed | perl ~/my_program3/src/goldfish_project/WGD2.triplet.pl $wd/input/$sp/triplet.
done
