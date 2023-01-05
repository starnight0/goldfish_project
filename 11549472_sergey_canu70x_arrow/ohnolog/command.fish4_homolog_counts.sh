############################################################
# count RBH genes
# CWD: fish4.cluster/czl_ohno_syn.out3
############################################################
cat rescue_m1.6.cluster.txt | awk -F$'\t' '$3~/ZF/ && $3~/GF:2/ {print $4}' | perl -ne '
chomp; my @t=split /;/;
my %g;
for (my $i=0;$i<@t; $i++) {
    foreach my $aa (split /,/, $t[$i]) {
        my @aa1;
        if ($aa=~m/^\((.*)\)$/) { @aa1 = split /\|/, $1; } else { @aa1=($aa); }
        my @gs;
        my $sp0;
        foreach my $aa1 (@aa1) {
            my ($sp,$chr,$id,@name) = split /:/, $aa1;
            my $name = join(":",@name);
            push @gs, [$chr,$id,$name];
            if (!defined $sp0) { $sp0 = $sp; }
        }
        push @{$g{$sp0}}, \@gs;
    }
}
for (my $i=1; $i<@{$g{GF}}; $i++) {
    my $g0 = $g{GF}[0];
    my $g1 = $g{GF}[$i];
    for (my $i0=0; $i0<@$g0; $i0++) {
        my $gg0 = $g0->[$i0];
        for (my $i1=0; $i1<@$g1; $i1++) {
            my $gg1 = $g1->[$i1];
            print join("\t", @{$gg0}), "\t", join("\t", @{$gg1}), "\t$i0\t$i1\n";
        }
    }
} ' > rescue_m1.6.cluster.pair.0.txt   
cat rescue_m1.6.cluster.txt | awk -F$'\t' '$3~/ZF/ && $3~/GF:2/ {print $4}' | perl -ne '
chomp; my @t=split /;/;
my %g;
for (my $i=0;$i<@t; $i++) {
    foreach my $aa (split /,/, $t[$i]) {
        my @aa1;
        if ($aa=~m/^\((.*)\)$/) { @aa1 = split /\|/, $1; } else { @aa1=($aa); }
        my @gs;
        my $sp0;
        foreach my $aa1 (@aa1) {
            my ($sp,$chr,$id,@name) = split /:/, $aa1;
            my $name = join(":",@name);
            push @gs, [$chr,$id,$name];
            if (!defined $sp0) { $sp0 = $sp; }
        }
        push @{$g{$sp0}}, \@gs;
    }
}
my $zf_chr = $g{ZF}[0][0][0];
my $gf_chr1 = "LG".$zf_chr;
my $gf_chr2 = "LG".($zf_chr+25);
for (my $i=1; $i<@{$g{GF}}; $i++) {
    my $g0 = $g{GF}[0];
    my $g1 = $g{GF}[$i];
    for (my $i0=0; $i0<@$g0; $i0++) {
        my $gg0 = $g0->[$i0];
        for (my $i1=0; $i1<@$g1; $i1++) {
            my $gg1 = $g1->[$i1];
            if ($gg0->[0] eq $gf_chr1 && $gg1->[0] eq $gf_chr2 || $gg0->[0] eq $gf_chr2 && $gg1->[0] eq $gf_chr1) { 
                print join("\t", @{$gg0}), "\t", join("\t", @{$gg1}), "\t$i0\t$i1\n";
            }
        }
    }
} ' > rescue_m1.6.cluster.pair.txt   

cat rescue_m1.6.cluster.pair.txt | awk '{print $2"\t"$0}' | sort -k1,1 > a;
join -t$'\t' /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.gene.unmasked.gids a  | cut -f 2- |  awk '{print $5"\t"$0}' | sort -k1,1 > a1;
join -t$'\t' /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.gene.unmasked.gids a1 | cut -f 2- | sort -k1,1 | sort -k2,2 > rescue_m1.6.cluster.pair.norm_coverage.txt   ;

zcat ../../fish17.blastn/pairs.gene.f3/GF.GF.f3.join.m6.gz | perl -ne '
BEGIN {
    open IN, "<rescue_m1.6.cluster.pair.norm_coverage.txt" or die; 
    while(<IN>) {
        if (m/^#/ || m/^\s*$/) {next;} chomp;
        my @t = split /\t/; $pair{$t[1]}{$t[4]}=\@t;
    }
    close IN;
}
chomp; my @t=split /\t/;
if ($t[0]=~m/^GF\|(.*)$/) { $t[0] = $1; }
if ($t[1]=~m/^GF\|(.*)$/) { $t[1] = $1; }
if (exists $pair{$t[0]}{$t[1]} && !exists $pass{$t[0]}{$t[1]}) {
    $pass{$t[0]}{$t[1]} = 1;
    my $p = $pair{$t[0]}{$t[1]};
    $p->[6]=$t[2];
    $p->[7]=$t[11];
    $p->[8]=$t[18];
    $p->[9]=$t[19];
    if ($t[18]>=50 && $t[19]>=50) { print join("\t", @$p), "\n"; }
} elsif (exists $pair{$t[1]}{$t[0]} && !exists $pass{$t[1]}{$t[0]}) {
    $pass{$t[1]}{$t[0]} = 1;
    my $p = $pair{$t[1]}{$t[0]};
    $p->[6]=$t[2];
    $p->[7]=$t[11];
    $p->[8]=$t[18];
    $p->[9]=$t[19];
    if ($t[18]>=50 && $t[19]>=50) { print join("\t", @$p), "\n"; }
} ' > rescue_m1.6.cluster.pair_extend.norm_coverage.txt

###########################################
cp ../../../carAur03/big/WGD/gene_pairs.1.txt ./pair_from_chainnet.txt
zcat ../../fish17.blastn/pairs.gene.f3/GF.GF.f3.join.m6.gz | perl -ne '
BEGIN {
    open IN, "<pair_from_chainnet.txt" or die; 
    while(<IN>) {
        if (m/^#/ || m/^\s*$/) {next;} chomp;
        my @t = split /\t/; $pair{$t[1]}{$t[4]}=\@t;
    }
    close IN;
}
chomp; my @t=split /\t/;
if ($t[0]=~m/^GF\|(.*)$/) { $t[0] = $1; }
if ($t[1]=~m/^GF\|(.*)$/) { $t[1] = $1; }
if (exists $pair{$t[0]}{$t[1]} && !exists $pass{$t[0]}{$t[1]}) {
    $pass{$t[0]}{$t[1]} = 1;
    my $p = $pair{$t[0]}{$t[1]};
    $p->[6]=$t[2];
    $p->[7]=$t[11];
    $p->[8]=$t[18];
    $p->[9]=$t[19];
    if ($t[18]>=50 && $t[19]>=50) { print join("\t", @$p), "\n"; }
} elsif (exists $pair{$t[1]}{$t[0]} && !exists $pass{$t[1]}{$t[0]}) {
    $pass{$t[1]}{$t[0]} = 1;
    my $p = $pair{$t[1]}{$t[0]};
    $p->[6]=$t[2];
    $p->[7]=$t[11];
    $p->[8]=$t[18];
    $p->[9]=$t[19];
    if ($t[18]>=50 && $t[19]>=50) { print join("\t", @$p), "\n"; }
} ' > pair_from_chainnet.extended.txt

######################################################


mkdir stat
cd stat
sps="ZF GC CC GF"
for i1 in `seq 1 4`
do
    sp1=`echo $sps | cut -d' ' -f $i1`; j1=`expr $i1 + 1`;
for i2 in `seq $j1 4`
do
    sp2=`echo $sps | cut -d' ' -f $i2`; j2=`expr $i2 + 1`;
    in=../RBH.rescue.pairwise.$sp1.$sp2.anchor_to_cluster.txt
    for sp in $sp1 $sp2
    do
        cat $in | grep '^'$sp | cut -f 3 | sort -k1,1 > tmp.${sp1}_$sp2.$sp
    done
done
done
cat tmp.ZF_GC.ZF tmp.ZF_CC.ZF | sort | uniq -d > tmp.ZF_GC_CC.ZF
cat tmp.ZF_GC.ZF tmp.ZF_GF.ZF | sort | uniq -d > tmp.ZF_GC_GF.ZF
cat tmp.ZF_GF.ZF tmp.ZF_CC.ZF | sort | uniq -d > tmp.ZF_CC_GF.ZF

cat tmp.ZF_GC.GC tmp.GC_CC.GC | sort | uniq -d > tmp.ZF_GC_CC.GC
cat tmp.ZF_GC.GC tmp.GC_GF.GC | sort | uniq -d > tmp.ZF_GC_GF.GC
cat tmp.GC_CC.GC tmp.GC_GF.GC | sort | uniq -d > tmp.GC_CC_GF.GC

cat tmp.ZF_CC.CC tmp.GC_CC.CC | sort | uniq -d > tmp.ZF_GC_CC.CC
cat tmp.ZF_CC.CC tmp.CC_GF.CC | sort | uniq -d > tmp.ZF_CC_GF.CC
cat tmp.GC_CC.CC tmp.CC_GF.CC | sort | uniq -d > tmp.GC_CC_GF.CC

cat tmp.ZF_GF.GF tmp.GC_GF.GF | sort | uniq -d > tmp.ZF_GC_GF.GF
cat tmp.ZF_GF.GF tmp.CC_GF.GF | sort | uniq -d > tmp.ZF_CC_GF.GF
cat tmp.GC_GF.GF tmp.CC_GF.GF | sort | uniq -d > tmp.GC_CC_GF.GF

cat tmp.??_??_??.ZF | sort | uniq > tmp.3.ZF
cat tmp.??_??_??.GC | sort | uniq > tmp.3.GC
cat tmp.??_??_??.CC | sort | uniq > tmp.3.CC
cat tmp.??_??_??.GF | sort | uniq > tmp.3.GF

cat tmp.ZF_GC_CC.ZF tmp.ZF_GC_GF.ZF | sort | uniq -d > ZF_GC_CC_GF.ZF
cat tmp.ZF_GC_CC.GC tmp.ZF_GC_GF.GC | sort | uniq -d > ZF_GC_CC_GF.GC
cat tmp.ZF_GC_CC.CC tmp.GC_CC_GF.CC | sort | uniq -d > ZF_GC_CC_GF.CC
cat tmp.ZF_GC_GF.GF tmp.ZF_CC_GF.GF | sort | uniq -d > ZF_GC_CC_GF.GF

for sps in ZF_GC_CC ZF_GC_GF ZF_CC_GF GC_CC_GF
do
    for sp in `echo $sps | sed 's/_/ /g'`
    do
        cat tmp.$sps.$sp ZF_GC_CC_GF.$sp | sort | uniq -d > a
        cat tmp.$sps.$sp a | sort | uniq -u > $sps.$sp
    done
done

for sps in ZF_GC ZF_CC ZF_GF GC_CC GC_GF CC_GF
do
    for sp in `echo $sps | sed 's/_/ /g'`
    do
        cat tmp.$sps.$sp tmp.3.$sp | sort | uniq -d > a;
        cat tmp.$sps.$sp a | sort | uniq -u > $sps.$sp;
    done
done

wc -l [A-Z]*.?? | awk '{print $2"\t"$1}' | grep -v total | perl -ne 'chomp; my @t=split /\t/; 
my ($sps,$sp)=split /\./, $t[0];
if (!exists $count{$sps}) { $count{$sps}={ZF=>0,GC=>0,CC=>0,GF=>0}; }
$count{$sps}{$sp}+=$t[1];
END {
    print join "\t", qw(SPS ZF GC CC GF), "\n"; 
    foreach my $sps (sort keys(%count)) {print join "\t", ($sps,$count{$sps}{ZF},$count{$sps}{GC},$count{$sps}{CC},$count{$sps}{GF}), "\n"}}' > RBH.stat.txt

###############################################
# homolog stat
###############################################
cd fish17.blastn
mkdir stat
cd stat
sps="ZF GC CC GF"
for i1 in `seq 1 4`
do
    sp1=`echo $sps | cut -d' ' -f $i1`; j1=`expr $i1 + 1`;
for i2 in `seq $j1 4`
do
    sp2=`echo $sps | cut -d' ' -f $i2`; j2=`expr $i2 + 1`;
    in=../pairs.gene.f3/$sp1.$sp2.f3.join.m6.gz
    zcat $in | cut -f 1 | cut -d'|' -f 2 | sort -k1,1 | uniq > tmp.${sp1}_$sp2.$sp1
    zcat $in | cut -f 2 | cut -d'|' -f 2 | sort -k1,1 | uniq > tmp.${sp1}_$sp2.$sp2
done
done
cat tmp.ZF_GC.ZF tmp.ZF_CC.ZF | sort | uniq -d > tmp.ZF_GC_CC.ZF
cat tmp.ZF_GC.ZF tmp.ZF_GF.ZF | sort | uniq -d > tmp.ZF_GC_GF.ZF
cat tmp.ZF_GF.ZF tmp.ZF_CC.ZF | sort | uniq -d > tmp.ZF_CC_GF.ZF

cat tmp.ZF_GC.GC tmp.GC_CC.GC | sort | uniq -d > tmp.ZF_GC_CC.GC
cat tmp.ZF_GC.GC tmp.GC_GF.GC | sort | uniq -d > tmp.ZF_GC_GF.GC
cat tmp.GC_CC.GC tmp.GC_GF.GC | sort | uniq -d > tmp.GC_CC_GF.GC

cat tmp.ZF_CC.CC tmp.GC_CC.CC | sort | uniq -d > tmp.ZF_GC_CC.CC
cat tmp.ZF_CC.CC tmp.CC_GF.CC | sort | uniq -d > tmp.ZF_CC_GF.CC
cat tmp.GC_CC.CC tmp.CC_GF.CC | sort | uniq -d > tmp.GC_CC_GF.CC

cat tmp.ZF_GF.GF tmp.GC_GF.GF | sort | uniq -d > tmp.ZF_GC_GF.GF
cat tmp.ZF_GF.GF tmp.CC_GF.GF | sort | uniq -d > tmp.ZF_CC_GF.GF
cat tmp.GC_GF.GF tmp.CC_GF.GF | sort | uniq -d > tmp.GC_CC_GF.GF

cat tmp.??_??_??.ZF | sort | uniq > tmp.3.ZF
cat tmp.??_??_??.GC | sort | uniq > tmp.3.GC
cat tmp.??_??_??.CC | sort | uniq > tmp.3.CC
cat tmp.??_??_??.GF | sort | uniq > tmp.3.GF

cat tmp.ZF_GC_CC.ZF tmp.ZF_GC_GF.ZF | sort | uniq -d > ZF_GC_CC_GF.ZF
cat tmp.ZF_GC_CC.GC tmp.ZF_GC_GF.GC | sort | uniq -d > ZF_GC_CC_GF.GC
cat tmp.ZF_GC_CC.CC tmp.GC_CC_GF.CC | sort | uniq -d > ZF_GC_CC_GF.CC
cat tmp.ZF_GC_GF.GF tmp.ZF_CC_GF.GF | sort | uniq -d > ZF_GC_CC_GF.GF

for sps in ZF_GC_CC ZF_GC_GF ZF_CC_GF GC_CC_GF
do
    for sp in `echo $sps | sed 's/_/ /g'`
    do
        cat tmp.$sps.$sp ZF_GC_CC_GF.$sp | sort | uniq -d > a
        cat tmp.$sps.$sp a | sort | uniq -u > $sps.$sp
    done
done

for sps in ZF_GC ZF_CC ZF_GF GC_CC GC_GF CC_GF
do
    for sp in `echo $sps | sed 's/_/ /g'`
    do
        cat tmp.$sps.$sp tmp.3.$sp | sort | uniq -d > a;
        cat tmp.$sps.$sp a | sort | uniq -u > $sps.$sp;
    done
done

wc -l [A-Z]*.?? | awk '{print $2"\t"$1}' | grep -v total | perl -ne 'chomp; my @t=split /\t/; 
my ($sps,$sp)=split /\./, $t[0];
if (!exists $count{$sps}) { $count{$sps}={ZF=>0,GC=>0,CC=>0,GF=>0}; }
$count{$sps}{$sp}+=$t[1];
END {
    print join "\t", qw(SPS ZF GC CC GF), "\n"; 
    foreach my $sps (sort keys(%count)) {print join "\t", ($sps,$count{$sps}{ZF},$count{$sps}{GC},$count{$sps}{CC},$count{$sps}{GF}), "\n"}}' > homo_gene.f3.stat.txt

###################
# other
###################
cat rescue_m1.6.cluster.txt | perl -e '
open IN, "zcat ../../fish17.blastn/pairs.gene.f3/ZF.GF.f3.join.m6.gz |";
my %pair;
while(<IN>) {
    chomp; my @t=split /\t/;
    my @t1=split /\|/,$t[0];   my @t2=split /\|/,$t[1];
    $pair{$t1[1]}{$t2[1]} = $t[2];
}
close IN;
while(<>) {
    my %g;
    chomp; my @t=split /\t/;
    if ($t[2]=~m/ZF:/ && $t[2]=~m/GF:/) {
        my @t1 = split /;/, $t[3];   my %gids;
        foreach my $t1 (@t1) {
            foreach my $t2 (split /,/, $t1) {
    my $aa = $t2;   my @aa1;   my @gs;   my $sp0;
    if ($aa=~m/^\((.*)\)$/) { @aa1 = split /\|/, $1; } else { @aa1=($aa); }
    foreach my $aa1 (@aa1) {
        my ($sp,$chr,$id,@name) = split /:/, $aa1;
        my $name = join(":",@name);
        push @gs, [$chr,$id,$name];
        if (!defined $sp0) { $sp0 = $sp; }
    }
    push @{$g{$sp0}}, \@gs;
            }
        }
        my $sp1="ZF"; my $sp2="GF";
        for (my $i=0; $i<@{$g{ZF}}; $i++) {
        for (my $j=0; $j<@{$g{GF}}; $j++) {
            my $g0 = $g{ZF}[$i];   my $g1 = $g{GF}[$j];
            for (my $i0=0; $i0<@$g0; $i0++) {
                my $gg0 = $g0->[$i0];
                for (my $i1=0; $i1<@$g1; $i1++) {
                    my $gg1 = $g1->[$i1];
                    my $gid0 = $gg0->[1];
                    my $gid1 = $gg1->[1];
                    if (exists $pair{$gid0}{$gid1}) {
                        print $gid0, "\t", $gid1, "\t", $pair{$gid0}{$gid1}, "\n";
                    } else {
                        print $gid0, "\t", $gid1, "\t", 0, "\n";
                    }
                }
            }
        } }
    }
}
' > tmp.ZF.GF
 
