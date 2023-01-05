#!/bin/sh
## sbatch --gres=lscratch:200 --time=168:00:00 --mem=16g -c 8 -a 1-684 --job-name run04.maker -o log2/run04.maker.o%A.%a run04.maker.sh
# biowulf module
export TMPDIR=/data/chenz11/tmp/maker4a
module load bamtools exonerate trf
# module repeatmasker 
export PATH=$PATH:$HOME/program/snap:$HOME/program/maker1/bin:$HOME/program/augustus-3.2.2/bin:$HOME/program/augustus-3.2.2/scripts:$HOME/program/EVidenceModeler-1.1.1:$HOME/program/EVidenceModeler-1.1.1/EvmUtils:/home/chenz11/program/gm_et_linux_64/gmes_petap:$HOME/program/RAPSearch2.24_64bits/bin:$HOME/program/RepeatMasker:$HOME/bin
export PERL5LIB=$HOME/lib/perl5/site_perl/5.18.2:$HOME/program/maker1/perl/lib:$HOME/program/maker1/lib:$HOME/program/EVidenceModeler-1.1.1/PerlLib

i=$SLURM_ARRAY_TASK_ID
#genome=`head -n $i ../by_2Mbp.list | tail -n 1` 
genome=../by_ctg/tig00004487_arrow.fa
maker -c $SLURM_CPUS_PER_TASK -base goldfish.arrow.renamed  -g $genome
