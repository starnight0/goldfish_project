cwd=`pwd`
datadir=~/data/goldfish/11549472/sergey_canu70x/arrow/miRBase_run1
genome=~/data/goldfish/11549472/sergey_canu70x/arrow/goldfish.arrow.renamed.masked.fasta
genome_name=goldfish.arrow.renamed.masked
file=$datadir/r01.blastn.miRNA.sh
mir=~/data/miRBase/hairpin.fa
mir1=~/data/miRBase/hairpin.U2T.fa
echo '#!/bin/sh' > $file
echo 'module load blast/2.6.0+ exonerate' >> $file
echo "if ! [ -f goldfish.arrow.renamed.masked.miRNA.m6 ]; then blastn -query $mir -db $genome -evalue 0.001 -perc_identity 90 -word_size 9 -num_threads \$SLURM_CPUS_PER_TASK -outfmt \"6 std btop qlen slen score gaps\" -num_alignments 100000 -out goldfish.arrow.renamed.masked.miRNA.m6; fi" >> $file
echo "if ! [ -f goldfish.arrow.renamed.masked.miRNA.m6.region ]; then perl ~/my_program3/src/assembly/fetch_blastn_region_for_exonerate.pl -i goldfish.arrow.renamed.masked.miRNA.m6 -o goldfish.arrow.renamed.masked.miRNA.m6.region; fi" >> $file 

odir=`pwd`
cd $datadir
if ! [ -d ex_split ]
then 
mkdir ex_split
split -a 3 -l 1000 -d $datadir/goldfish.arrow.renamed.masked.miRNA.m6.region ex_split/in.
file=$datadir/r02.miRNA.ex.sh
ls ex_split/in.* | sed 's/^.*\.\([0-9]\+\)$/\1/' > ex_split.input.i
fi
cd $odir
echo '#!/bin/sh' > $file
echo '#SBATCH --time=4:00:00' >> $file
echo '#SBATCH -p quick' >> $file
echo '#SBATCH -o r02.miRNA.ex.o%A.%a' >> $file
echo '#SBATCH -J r02.miRNA.ex' >> $file
echo 'module load blast/2.6.0+ exonerate' >> $file
echo 'i=$SLURM_ARRAY_TASK_ID' >> $file
echo 'j=`head -n $i ex_split.input.i | tail -n 1`' >> $file
echo 'perl ~/my_program3/src/assembly/czl_exonerate_after_blasn.pl -i ex_split/in.$j -o ex_split/out.$j -q '$mir1' -t '$genome' -p $SLURM_CPUS_PER_TASK -fl 1000' >> $file
# echo "perl ~/my_program3/src/assembly/czl_exonerate_after_blasn.pl -i goldfish.arrow.renamed.masked.miRNA.m6.region -o goldfish.arrow.renamed.masked.miRNA.ex -q $mir1 -t $genome -p \$SLURM_CPUS_PER_TASK -fl 1000" >> $file

file=$datadir/r01.rfam.sh
dir=$datadir/rfam_by_2Mbp
if ! [ -d $dir ]; then mkdir $dir; fi
echo '#!/bin/sh' > $file
echo '#SBATCH --time=4:00:00' >> $file
echo '#SBATCH -p quick' >> $file
echo '#SBATCH -o rfam_by_2Mbp/r01.rfam.o%A.%a' >> $file
echo '#SBATCH -J r01.rfam' >> $file
echo 'i=$SLURM_ARRAY_TASK_ID' >> $file
echo 'in=`cat ../by_2Mbp.list | head -n $i | tail -n 1`' >> $file
echo 'out=rfam_by_2Mbp/$i' >>  $file
echo '~/program/infernal-1.1.2/bin/cmscan --rfam --cut_ga --nohmmonly --tblout $out.cmscan.tblout --fmt 2 --clanin ~/data/RFAM/12.3/Rfam.clanin --cpu $SLURM_CPUS_PER_TASK ~/data/RFAM/12.3/Rfam.1_1.cm $in > $out.cmscan' >> $file

file=$datadir/r01.tRNAscan.sh
echo '#!/bin/sh' > $file
echo 'export PATH=$PATH:$HOME/bin' >> $file
echo 'export PERL5LIB=$PERL5LIB:$HOME/bin' >> $file
echo "tRNAscan-SE -o $genome_name.tRNA.out -f $genome_name.tRNA.structure -m $genome_name.tRNA.summary $genome" >>$file

cd $datadir
file=r02.miRNA.rnafold.sh
echo '#!/bin/sh' > $file
echo '#SBATCH --time=4:00:00' >> $file
echo '#SBATCH --mem=8g' >> $file
echo '#SBATCH -p quick' >> $file
echo '#SBATCH -o ex_split/r02.miRNA.rnafold.o%A' >> $file
echo '#SBATCH -J r02.miRNA.rnafold' >> $file
echo 'module load viennarna' >> $file
echo 'RNAfold -p -d2 --noLP --noPS < ex_split/all.out.fasta > ex_split/all.out.rnafold' >> $file
cd $cwd

# run miRNA.sh

