##############################
# remove heterozygous contigs for carAur03
##############################
# CWD: carAur03.vs.carAur03
cd carAur03.vs.carAur03
cat all.target.syn.net.chain | perl -e '
my %mask;
open IN, "<../../carAur01/carAur01.masked_ctg_id";
while(<IN>) { if (m/^#/ || m/^\s*$/) {next;} chomp; $mask{$_}++; }
close IN;
my @a;
my $filt=0;
while(<>) {
    if (m/^#/) {next;} chomp;
    if (m/^chain/) {
        my @t = split / /, $_;
        if (exists $mask{$t[2]} || exists $mask{$t[7]}) { $filt=1; @a=();}
        else { $filt=0; @a=($_); }
    } elsif (m/\S/) {
        if (!$filt) { push @a, $_; }
    } else {
        if (!$filt) {
            foreach my $a (@a) { print $a, "\n"; }
            print "\n";
        }
    }
}' > all.target.syn.net.no_het.chain

chainToAxt -maxGap=20 all.target.syn.net.no_het.chain ../../carAur03/carAur03.2bit  ../../carAur03/carAur03.2bit stdout | axtFilter -minScore=1000 stdin | axtToChain stdin ../../carAur03/carAur03.sm.fa.fai ../../carAur03/carAur03.sm.fa.fai stdout | grep '^chain' | awk -v OFS=$'\t' -v n=0 '{ if ($5=="+") {print $3,$6,$7,n;} else {print $3,$4-$7,$4-$6,n} if ($10=="+") {print $8,$11,$12,n} else {print $8,$9-$12,$9-$11,n} n++}' > all.target.syn.net.no_het.chain.20.bed4
#chainSwap all.target.syn.net.no_het.chain stdout | chainToAxt -maxGap=20 stdin ../../carAur03/carAur03.2bit  ../../carAur03/carAur03.2bit stdout | axtFilter -minScore=1000 stdin | axtToChain stdin ../../carAur03/carAur03.sm.fa.fai ../../carAur03/carAur03.sm.fa.fai stdout | grep '^chain' | awk -v OFS=$'\t' -v n=0 '{ if ($5=="+") {print $3,$6,$7,n;} else {print $3,$4-$7,$4-$6,n} if ($10=="+") {print $8,$11,$12,n} else {print $8,$9-$12,$9-$11,n} n++}' >> all.target.syn.net.no_het.chain.20.bed4
cat all.target.syn.net.no_het.chain.20.bed4 | cut -f 1-3 | sort -k1,1 -k2,2n > all.target.syn.net.no_het.chain.20.bed3
bedtools merge -i all.target.syn.net.no_het.chain.20.bed3 | sort -k1,1 -k2,2n > all.target.syn.net.no_het.chain.20.merged.bed3

# CWD: carp_ncbi.vs.carp_ncbi
cd carp_ncbi.vs.carp_ncbi
chainToAxt -maxGap=20 all.target.syn.net.chain ~/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.2bit  ~/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.2bit stdout | axtFilter -minScore=1000 stdin | axtToChain stdin  ~/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.sm.fa.fai ~/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.sm.fa.fai stdout | grep '^chain' | awk -v OFS=$'\t' -v n=0 '{ if ($5=="+") {print $3,$6,$7,n;} else {print $3,$4-$7,$4-$6,n} if ($10=="+") {print $8,$11,$12,n} else {print $8,$9-$12,$9-$11,n} n++}' > all.target.syn.net.chain.20.bed4
#chainSwap all.target.syn.net.chain stdout | chainToAxt -maxGap=20 stdin ~/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.2bit  ~/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.2bit stdout | axtFilter -minScore=1000 stdin | axtToChain stdin  ~/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.sm.fa.fai ~/data/common_carp/NCBI_Cyprinus_carpio/Cyprinus_carpio.sm.fa.fai stdout | grep '^chain' | awk -v OFS=$'\t' -v n=0 '{ if ($5=="+") {print $3,$6,$7,n;} else {print $3,$4-$7,$4-$6,n} if ($10=="+") {print $8,$11,$12,n} else {print $8,$9-$12,$9-$11,n} n++}' >> all.target.syn.net.chain.20.bed4
cat all.target.syn.net.chain.20.bed4 | cut -f 1-3 | sort -k1,1 -k2,2n > all.target.syn.net.chain.20.bed3
bedtools merge -i all.target.syn.net.chain.20.bed3 | sort -k1,1 -k2,2n > all.target.syn.net.chain.20.merged.bed3


