chain=big/WGD/GF.ZF.net.long.non_Tovl.chain
# query exons with ID=chr:begin:end:strand:gid:tid:eid
qexon=
# target/reference exons with ID=chr:begin:end:strand:gid:tid:eid
texon=
out_dir=

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
		print join("\t", ($chr, $t[3], $t[4], "$chr:$t[3]:$t[4]:$t[6]:$gene_id:$transcript_id:$exon_id", 0, $t[6])), "\n";
	}
	' | sort -k1,1 -k2,2n > $out_dir/texon.bed;
texon=$out_dir/texon.bed;

cat big/carAur03.noM.exon.bed | perl -ne '
s/\s+$//;
my @t=split /\t/;
$eid = $t[3];
$tid = $eid; $tid=~s/:.*$//;
$gid = $tid; $gid=~s/_R[0-9]+$//;
print join("\t", (@t[0..2],"$t[0]:$t[1]:$t[2]:$t[5]:$gid:$eid",$t[4],$t[5])), "\n";
' | sort -k1,1 -k2,2n > $out_dir/qexon.bed;
qexon=$out_dir/qexon.bed;
crossmap bed $chain $qexon $out_dir/qexon.map.bed
cat $out_dir/qexon.map.bed | sort -k1,1 -k2,2n > $out_dir/qexon.map.sorted.bed 
bedtools intersect -nonamecheck -wao -f 0.5 -a $out_dir/qexon.map.bed -b $texon  > $out_dir/qexon.map.texon.bed

cat $out_dir/qexon.map.texon.bed | perl -e '
my %join_gene;
my ($gid20, @gid1s, %gid1s);
my @tot=(0) x 2;
my ($m1,$m2,$m3,$m4)=(0,0,0,0);
while(<>) {
	chomp;
	my @t=split /\t/;
	if ($t[3] eq "." || $t[9] eq ".") { next; }
	my @id1 = split /:/, $t[3];
	my @id2 = split /:/, $t[9];
	if (@id2<2) { next; }
	my $gid1 = $id1[4];
	my $gid2 = $id2[4];
	if (!defined $gid20) { $gid20=$gid2;}
	elsif (defined $gid20 && $gid20 ne $gid2) {
		$tot[1]++; $tot[0]+=@gid1s;
		if (@gid1s>1) { process($gid20); }
		@gid1s=(); %gid1s=(); $gid20=$gid2;
	}
	if (!exists $gid1s{$gid1}) { push @gid1s, $gid1; }
	$gid1s{$gid1}[0]++;
	$gid1s{$gid1}[1]+=$t[2]-$t[1];
}
$tot[1]++; $tot[0]+=@gid1s;
if (@gid1s>1) {  process($gid20); }
print STDERR "$tot[0]\t$tot[1]\t$m1\t$m2\t$m3\t$m4\n";

sub process {
	$jgid = join("+",sort @gid1s);
	if (!exists $join_gene{$jgid}) {
		my $gid20 = shift;
		my $n=0; my $l=0; my $n2 =0;
		foreach my $gid (@gid1s) { $n+=$gid1s{$gid}[0]; $l+=$gid1s{$gid}[1];}
		foreach my $gid (@gid1s) {
			my $f = sprintf("%.3f", $gid1s{$gid}[1]/$l);
			if ($f>0.1 || ($gid1s{$gid}[1]>=30 && $gid1s{$gid}[0]>2)) {$n2++;}
		}
		foreach my $gid (@gid1s) {
			my $f = sprintf("%.3f", $gid1s{$gid}[1]/$l);
			print "$gid\t$gid20\t$gid1s{$gid}[0]\t$gid1s{$gid}[1]\t$f\n";
			if ($f>0.1 || $gid1s{$gid}[0]>2) {$n2++;}
		}
		if ($n2>1) { $m3+=@gid1s; $m4++;}
		$join_gene{$jgid} = $gid20;
		$m1+=@gid1s;   $m2++;
	}
}
' > $out_dir/split_gene.txt




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
		if ($info{gene_biotype}=~m/protein_coding/ && defined $transcript_id) { print "$transcript_id\t$l\n"; }
		$transcript_id = $info{transcript_id};
		$l = 0;
	} elsif ($t[2] eq "exon") {
		$l+=$t[4]-$t[3];
	}
	' | sort -k1,1 > $out_dir/ttranscript_size.txt;


cat big/carAur03.noM.gene.unmasked.gff | perl -ne '
	if (m/^#/ || m/^\s*$/) { next; }
	s/\s+$//;
	my @t = split /\t/;
    my %info;
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
	} elsif ($t[2] eq "transcript" || $t[2]=~m/RNA/i) {
		if (defined $transcript_id) { print "$transcript_id\t$l\n"; }
		$transcript_id = $info{ID};
		$l = 0;
	} elsif ($t[2] eq "exon") {
		$l+=$t[4]-$t[3];
	}
	' | sort -k1,1 > $out_dir/qtranscript_size.txt;

