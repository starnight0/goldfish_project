datadir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/maker4a
cwd=`pwd`
cd $datadir

gunzip -c ~/data/ensembl85/ens85_fish_cc_gc.pep.short_name.fasta.gz > ens85_fish_cc_gc_gf.pep.fasta
cat carAur01.all.maker.proteins.short_name.fasta | sed '/^>/ s/^>/>carAur./'>> ens85_fish_cc_gc_gf.pep.fasta
gunzip -c ~/data/ensembl85/ens85_fish_cc_gc.rna.short_name.fasta.gz > ens85_fish_cc_gc_gf.rna.fasta
cat carAur01.all.maker.transcripts.short_name.fasta | sed '/^>/ s/^>/>carAur./' >> ens85_fish_cc_gc_gf.rna.fasta

cat ~/data/ncbi/refseq/vertebrate_mammalian/rna/*.rna.fna ~/data/ncbi/refseq/vertebrate_other/rna/*.rna.fna carAur01.all.maker.transcripts.short_name.fasta | sed '/^>/ s/\s.*$//' | sed '/^>/ s/[^0-9Z-z]\+$//' > refseq_gf.rna.fasta
makeblastdb -in refseq_gf.rna.fasta -dbtype nucl

cat ~/data/ncbi/refseq/vertebrate_mammalian/protein/*.protein.faa ~/data/ncbi/refseq/vertebrate_other/protein/*.protein.faa carAur01.all.maker.proteins.short_name.fasta | sed '/^>/ s/\s.*$//' | sed '/^>/ s/[^0-9Z-z]\+$//' > refseq_gf.pep.fasta
makeblastdb -in refseq_gf.pep.fasta -dbtype prot


mkdir -p t.blast/
file=t.blast/blastn.sh
echo '#!/bin/sh' > $file
echo '#SBATCH --mem=16g' >> $file
echo '#SBATCH --time=120:00:00' >> $file
echo '#SBATCH -c 32' >> $file
echo '#SBATCH -o blastn.o%A' >> $file
echo '#SBATCH -J blastn' >> $file
echo 'module load blast/2.6.0+' >> $file
echo 'if ! [ -f out.m6 ]; then blastn -task blastn -query ../ens85_fish_cc_gc_gf.rna.fasta -db ../ens85_fish_cc_gc_gf.rna.fasta -evalue 0.0001 -num_threads $SLURM_CPUS_PER_TASK -out out.m6 -outfmt "6 std btop qlen slen score gaps positive" -num_alignments 1000; fi' >> $file

file=t.blast/blastn.ref.sh
echo '#!/bin/sh' > $file
echo '#SBATCH --mem=20g' >> $file
echo '#SBATCH --time=240:00:00' >> $file
echo '#SBATCH -c 32' >> $file
echo '#SBATCH -o blastn.ref.o%A' >> $file
echo '#SBATCH -J blastn.ref' >> $file
echo 'module load blast/2.6.0+' >> $file
echo 'if ! [ -f out.refseq.m6.gz ]; then blastn -task blastn -query ../refseq_gf.rna.fasta -db ../refseq_gf.rna.fasta -evalue 0.0001 -num_threads $SLURM_CPUS_PER_TASK -outfmt "6 std btop qlen slen score gaps positive" -num_alignments 200 | gzip -c > out.refseq.m6.gz ; fi' >> $file


mkdir -p p.blast/ 
file=p.blast/blastp.sh
echo '#!/bin/sh' > $file
echo '#SBATCH --mem=16g' >> $file
echo '#SBATCH --time=120:00:00' >> $file
echo '#SBATCH -c 32' >> $file
echo '#SBATCH -o blastp.o%A' >> $file
echo '#SBATCH -J blastp' >> $file
echo 'module load blast/2.6.0+' >> $file
echo 'if ! [ -f out.m6 ]; then blastp -query ../ens85_fish_cc_gc_gf.pep.fasta -db ../ens85_fish_cc_gc_gf.pep.fasta -evalue 0.0001 -num_threads $SLURM_CPUS_PER_TASK -out out.m6 -outfmt "6 std btop qlen slen score gaps positive" -num_alignments 1000; fi' >> $file

file=p.blast/blastp.ref.sh
echo '#!/bin/sh' > $file
echo '#SBATCH --mem=32g' >> $file
echo '#SBATCH --time=240:00:00' >> $file
echo '#SBATCH -c 32' >> $file
echo '#SBATCH -o blastp.ref.o%A' >> $file
echo '#SBATCH -J blastp.ref' >> $file
echo 'module load blast/2.6.0+' >> $file
echo 'if ! [ -f out.refseq.m6.gz ]; then blastp -query ../refseq_gf.pep.fasta -db ../refseq_gf.pep.fasta -evalue 0.0001 -num_threads $SLURM_CPUS_PER_TASK -outfmt "6 std btop qlen slen score gaps positive" -num_alignments 1000 | gzip -c > out.refseq.m6.gz ; fi' >> $file


# after finishing t.blast/blastn.sh
# run:
if [ "x" != "x" ]
then
	for dir1 in t.blast q.blast
		cat $dir1/out.m6 | perl -e 'while(<>) { @t=split "\t"; 
			$q=$t[0]; $t=$t[1]; 
			if ($q=~m/^ENS/) {$qsp=substr $q,0,6;} elsif ($q=~m/^CI/) {$qsp="CI";} elsif ($q=~m/^CAFS/) {$qsp="CC";} else {$qsp="GF";}
			if ($t=~m/^ENS/) {$tsp=substr $t,0,6;} elsif ($t=~m/^CI/) {$tsp="CI";} elsif ($t=~m/^CAFS/) {$tsp="CC";} else {$tsp="GF";}
			my $sp1 = $qsp;
			my $sp2 = $tsp;
			if ( ($sp1 cmp $sp2) > 0 ) { $sp1=$tsp; $sp2=$qsp; }
			if (!exists $fh{$sp1}{$sp2}) { open $fh{$sp1}{$sp2}, "| gzip -c > '"$dir1"'/pairs/$sp1.$sp2.m6.gz"; $fh{$sp2}{$sp1}=$fh{$sp1}{$sp2}; }
			print {$fh{$sp1}{$sp2}} $_;
		}
		foreach my $sp1 (keys(%fh)) {
			foreach my $sp2 (keys(%{$fh{$sp1}})) {
				close $fh{$sp1}{$sp2};
			}
		}'
		gzip $dir1/out.m6
fi

cd $cwd


