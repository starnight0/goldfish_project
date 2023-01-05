# map ids
>r01.map_id.sh
echo '#!/bin/sh
#SBATCH -p quick --mem=8g --time=1:00:00 -o r01.map_id.o%A
~/program/maker/bin/maker_map_ids --prefix CARNA --suffix _R --iterate 1 carAur01.ncrna.gff > carAur01.ncrna.gff.maker_map_ids
cp carAur01.ncrna.gff carAur01.ncrna.old_id.gff; gzip  carAur01.ncrna.old_id.gff;
~/program/maker/bin/map_gff_ids carAur01.ncrna.gff.maker_map_ids carAur01.ncrna.gff
sed -i "s/Name=;//" carAur01.ncrna.gff
gffread  carAur01.ncrna.gff -g ../carAur01/carAur01.fa -w carAur01.ncrna.fa >stdout 2>stderr
' >> r01.map_id.sh



ens_gff=
gff1=../ncrna_run1/all.ex1.f.cuffcompare.combined.gtf
gff2=../miRBase_run1/all.ex.f.cuffcompare.combined.gtf.not_exon.gtf
gff3=../miRBase_run1/carAur01.rfam.cmscan.f.gff
gff4=../miRBase_run1/all.ex.out.f.gtf
out=all4.cuffcompare.combined.add_id.gff
cat ../ncrna_run1/ens85_ncrna_70.ex1/all.out.rmdup.gff | grep ENSDAR > ens85_ncrna_70.ex1.ZF.gff
cat ../ncrna_run1/ens85_ncrna_70.ex1/all.out.f.rmdup.gff | grep -v ENSDAR | sort -k1,1 -s -S 1G > ens85_ncrna_70.ex1.not_ZF.gff
cat ../ncrna_run1/ens85_ncrna_70.ex1/all.out.rmdup.gff | grep ENSDAR | sort -k1,1 -s -S 1G > ens85_ncrna_70.ex1.ZF.gff
cat ../miRBase_run1/carAur01.rfam.cmscan.f.gff | sort -k1,1 -s -S 1G > rfam.f.gff
cat ../miRBase_run1/all.ex.out.f.gtf | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > miR.f.gtf
cat ../ncrna_run1/NONCODE.ex1/all.out.f.rmdup.gff | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > NONCODE.f.gff
cat ../ncrna_run1/rnacentral.[1-9].ex1/all.out.f.rmdup.gff | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > rnacentral.1_9.f.gff
cat ../ncrna_run1/rnacentral.1[0-9].ex1/all.out.f.rmdup.gff | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > rnacentral.10_19.f.gff
cat ../ncrna_run1/rnacentral.2[0-9].ex1/all.out.f.rmdup.gff | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > rnacentral.20_29.f.gff
cat ../ncrna_run1/rnacentral.3[0-9].ex1/all.out.f.rmdup.gff | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > rnacentral.30_39.f.gff
cat ../ncrna_run1/rnacentral.4[0-9].ex1/all.out.f.rmdup.gff | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > rnacentral.40_49.f.gff
cat ../ncrna_run1/rnacentral.5[0-9].ex1/all.out.f.rmdup.gff | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > rnacentral.50_59.f.gff
cat ../ncrna_run1/rnacentral.6[0-8].ex1/all.out.f.rmdup.gff | awk '$3~/transcript/ || $3~/exon/' | sort -k1,1 -s -S 1G > rnacentral.60_68.f.gff


file=r02.add_name.sh
>$file
echo '#!/bin/sh' >> $file
echo 'cgff=all4.cuffcompare.combined.gtf
~/my_program3/src/utility/czl_gff_util.pl annot -i $cgff -i2 ens85_ncrna_70.ex1.ZF.gff -attr "ID=>zf_ens_id" -o all4.tmp1.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp1.gff -i2 ens85_ncrna_70.ex1.not_ZF.gff -attr "ID=>other_ens_id" -o all4.tmp2.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp2.gff -i2 rfam.f.gff -attr "Name=>rfam_name,rfam_id=>rfam_id,rfam_clan=>rfam_clan" -o all4.tmp3.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp3.gff -i2 miR.f.gtf -attr "ID=>mir_id" -o all4.tmp4.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp4.gff -i2 NONCODE.f.gff -attr "ID=>noncode_id" -o all4.tmp5.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp5.gff -i2 rnacentral.1_9.f.gff -attr "ID=>rnacentral_id" -o all4.tmp6.gff -ist;
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp6.gff -i2 rnacentral.10_19.f.gff -attr "ID=>rnacentral_id" -o all4.tmp7.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp7.gff -i2 rnacentral.20_29.f.gff -attr "ID=>rnacentral_id" -o all4.tmp8.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp8.gff -i2 rnacentral.30_39.f.gff -attr "ID=>rnacentral_id" -o all4.tmp9.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp9.gff -i2 rnacentral.40_49.f.gff -attr "ID=>rnacentral_id" -o all4.tmp10.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp10.gff -i2 rnacentral.50_59.f.gff -attr "ID=>rnacentral_id" -o all4.tmp11.gff -ist
~/my_program3/src/utility/czl_gff_util.pl annot -i all4.tmp11.gff -i2 rnacentral.60_68.f.gff -attr "ID=>rnacentral_id" -o all4.tmp12.gff -ist
' >> $file
cmd=r02.add_name; sbatch --mem=32G --time=4:00:00 -p quick -J $cmd -o $cmd.o%A $cmd.sh 
#cmd=r02.add_name; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m blast,TransDecoder,bedtools -g 16 --partition quick --time=4:00:00 -J $cmd --logdir=$cmd.log -m cufflinks -f $cmd.sh

cat all4.tmp5.gff | perl -e '
open IN, "</home/chenz11/data/ensembl85/all.85.gtpn";
while (<IN>) { 
if (m/^#/ || m/^\s*$/) { next; } chomp; my @t=split "\t"; 
if ($t[3] ne "." && $t[3]=~m/\S/) { $tr{$t[1]}=$t[3]; }
}
close IN;
while (<>) {
if (m/^#/ || m/^\s*$/) { print $_; next; } chomp; my @t=split "\t"; 
if ($t[2]=~m/gene|mRNA|transcript/) {
	my ($zf_ens_id, $ens_id, $mir_id, $rfam_name);
	foreach my $aa (split /\s*;\s*/, $t[8]) {
		my ($u,$v) = split /\s*=\s*/,$aa;
		$v=~s/^[0-9]+_//; $v=~s/^[0-9]+_//;
		if ($u=~m/zf_ens_id/) { $zf_ens_id=$v; }
		elsif ($u=~m/ens_id/) { $ens_id=$v; }
		elsif ($u=~m/^mir_id$/) { $mir_id=$v; }
		elsif ($u=~m/^rfam_name$/) { $rfam_name=$v; }
	}
	if (defined $zf_ens_id && exists $tr{$zf_ens_id}) {
		if ($t[8]=~m/;$/) { $t[8] .= "Name=$tr{$zf_ens_id}"; }
		else { $t[8] .= ";Name=$tr{$zf_ens_id}"; }
	}
	elsif (defined $ens_id && exists $tr{$ens_id}) {
		if ($t[8]=~m/;$/) { $t[8] .= "Name=$tr{$ens_id}"; }
		else { $t[8] .= ";Name=$tr{$ens_id}"; }
	}
	elsif (defined $mir_id && exists $tr{$mir_id}) {
		if ($t[8]=~m/;$/) { $t[8] .= "Name=$tr{$mir_id}"; }
		else { $t[8] .= ";Name=$tr{$mir_id}"; }
	}
	elsif (defined $rfam_name && exists $tr{$rfam_name}) {
		if ($t[8]=~m/;$/) { $t[8] .= "Name=$tr{$rfam_name}"; }
		else { $t[8] .= ";Name=$tr{$rfam_name}"; }
	}
}
	print join("\t",@t), "\n";
}
' > carAur01.ncrna.gff
