#!/bin/sh
module load bwa
bwa index carAur01.all.maker.transcripts.short_name.fasta
bwa index carAur01.all.evm.transcripts.short_name.fasta
