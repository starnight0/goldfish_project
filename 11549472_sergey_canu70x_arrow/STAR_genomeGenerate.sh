#!/bin/sh
#SBATCH -p quick --time=4:00:00 -c 4 -J STAR_genomeGenerate -o STAR_genomeGenerate.o%A --mem=64g
module load STAR
ref=carAur01.fasta
genome_dir=`pwd`/STAR_carAur01_maker
gff=maker4a/carAur01.gene.maker.gff
# for gff use --sjdbGTFtagExonParentTranscript Parent
STAR --runMode genomeGenerate \
	--runThreadN 4 \
	--genomeDir $genome_dir \
	--genomeFastaFiles $ref \
	--sjdbGTFfile $gff \
	--sjdbGTFtagExonParentTranscript Parent \
	--sjdbOverhang 120
gtf=maker4a/carAur01.gene.evm.gff.gz
genome_dir=`pwd`/STAR_carAur01_evm
STAR --runMode genomeGenerate \
	--runThreadN 4 \
	--genomeDir $genome_dir \
	--genomeFastaFiles $ref \
	--sjdbGTFfile $gff \
	--sjdbGTFtagExonParentTranscript Parent \
	--sjdbOverhang 120
