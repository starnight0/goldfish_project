
file=r01.split_fasta.sh
echo "Produce $file"
>$file
for sm in `cat sample`
do
    echo 'cd '$sm'.trinity; d=Trinity-GG.short_name; if ! [ -d $d.by_2Mbp ]; then mkdir $d.by_2Mbp; ~/my_program3/src/utility/czl_fasta_split -m 2 -n 2000000 -i $d.fasta -o $d.by_2Mbp/ ; fi' >> $file
done

file=r01.TransDecoder.sh
echo "Produce $file"
>$file
for sm in `cat sample`
do
    echo 'cd '$sm'.trinity; TransDecoder.LongOrfs  -t Trinity-GG.fasta -m 30 -S; cd ..' >> $file;
    echo 'cd '$sm'.trinity/Trinity-GG.fasta.transdecoder_dir; cat longest_orfs.pep | sed -e "/^>/ s/\s.*$//" -e "s/\*$//" > longest_orfs.no_stop_codon.pep; cd ../../' >> $file;
done
#cmd=r01.TransDecoder; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m blast,TransDecoder -g 4 --time=48:00:00 -J $cmd --logdir=$cmd.log -f $cmd.sh

file=r02.interpro.sh
echo "Produce $file"
>$file
for sm in `cat sample`
do
    echo 'cd '$sm'.trinity/Trinity-GG.fasta.transdecoder_dir; interproscan -t p -f TSV,GFF3,HTML -iprlookup --goterms --pathways --minsize 30 longest_orfs.no_stop_codon.pep interProScan5_run1 1000;' >> $file
done
#cmd=r02.interpro; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m blast,interproscan -g 4 --time=48:00:00 -J $cmd --logdir=$cmd.log -f $cmd.sh

IFS=$'\n'
for db in `cat prot_db`
do
    dbname=`echo $db | cut -f 1`;
    dbfile=`echo $db | cut -f 2`;
    file=r02.blastp.$dbname.sh
    echo "Produce $file"
    >$file
    for sm in `cat sample`
    do
        dir=$sm.trinity/Trinity-GG.fasta.transdecoder_dir;
        f=$dir/longest_orfs.no_stop_codon.pep
        echo 'blastp -db '$dbfile' -query '$f' -outfmt "6 std score qlen slen gaps positive" -evalue 1e-4 -num_threads 32 -max_target_seqs 100 | gzip -c > '$dir'/longest_orfs.pep.'$dbname'.blastp.m6.gz' >> $file;
    done
done
#cmd=r01.blastp; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm --gres=lscratch:50 -m blast -t 32 -g 56 --time=48:00:00 -J $cmd --logdir=$cmd.log -f $cmd.sh

IFS=$'\n'
for db in `cat nucl_db`
do
    dbname=`echo $db | cut -f 1`;
    dbfile=`echo $db | cut -f 2`;
    opt=`echo $db | cut -f 3`;
    file=r02.blastn.$dbname.sh
    echo "Produce $file"
    >$file
    for sm in `cat sample`
    do
        dir=$sm.trinity;
        f=$dir/Trinity-GG.short_name.fasta
        echo 'blastn -db '$dbfile' -query '$f' -outfmt "6 std score qlen slen gaps positive" -evalue 1e-4 -num_threads $SLURM_CPUS_PER_TASK '$opt' | gzip -c > '$dir'/Trinity-GG.'$dbname'.blastn.m6.gz' >> $file;
    done
done
#cmd=r01.blastn; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm --gres=lscratch:16 -t 32 -m blast -g 16 --time=96:00:00 -J $cmd --logdir=$cmd.log -f $cmd.sh


for sm in `cat sample`
do
    dir=$sm.trinity;
    cd $dir
    if ! [ -d blastx_out ]; then mkdir blastx_out; fi
    IFS=$'\n'
    file1=r02.blastx.cat.sh
    >$file1
    for db in `cat ../prot_db`
    do
        dbname=`echo $db | cut -f 1`;
        dbfile=`echo $db | cut -f 2`;
        file=r02.blastx.$dbname.sh
        echo "Produce $file"
        >$file
        if ! [ -d Trinity-GG.short_name.by_2Mbp ]; then break; fi
        for f in `ls -v Trinity-GG.short_name.by_2Mbp/*.fa`
        do
            i=`echo $f | sed -e 's/^.*\///' -e 's/\.fa$//'`
            if [ "$dbname" == "nr" ]
            then
                echo 'blastx -task blastx-fast -db '$dbfile' -query '$f' -outfmt "6 std score qlen slen gaps positive" -evalue 1e-4 -num_threads $SLURM_CPUS_PER_TASK -max_target_seqs 100 | gzip -c > blastx_out/'$i'.'$dbname'.blastx.m6.gz' >> $file;
            else
                echo 'blastx -db '$dbfile' -query '$f' -outfmt "6 std score qlen slen gaps positive" -evalue 1e-4 -num_threads $SLURM_CPUS_PER_TASK -max_target_seqs 100 | gzip -c > blastx_out/'$i'.'$dbname'.blastx.m6.gz' >> $file;
            fi
        done
        echo 'if ! [ -f Trinity-GG.'$dbname'.blastx.m6.gz ]; then  zcat blastx_out/*.'$dbname'.blastx.m6.gz | gzip -c > Trinity-GG.'$dbname'.blastx.m6.gz;  fi' >> $file1
    done
    cd ..
done
#cmd=r02.blastx; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -t 32 -m blast -g 24 --time=96:00:00 --merge-output -J $cmd --logdir=$cmd.log -f $cmd.sh

# run exonerate
file=r03.fetch_ex_region.sh
>$file
for sm in `cat sample`
do
    dir=$sm.trinity;
    exdir=$dir/Trinity-GG.carAur01.ex
    if ! [ -d $exdir ]; then mkdir $exdir; fi
    echo '~/my_program3/src/annot_genome/czl_fetch_blastn_region_for_exonerate.pl -i '$dir'/Trinity-GG.carAur01.blastn.m6.gz -o '$dir'/Trinity-GG.carAur01.blastn.region --less-qovl;    split -d -a 3 -l 500000 '$dir'/Trinity-GG.carAur01.blastn.region '$exdir'/' >> $file
done
#cmd=r03.fetch_ex_region; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m blast -g 32 --time=96:00:00 --merge-output -J $cmd --logdir=$cmd.log -f $cmd.sh

file1=r03.ex2.cat.sh
>$file1
for sm in `cat sample`
do
    dir=$sm.trinity;
    exdir=$dir/Trinity-GG.carAur01.ex2
    cd $exdir
    file=r03.ex.sh
    echo "Produce $file"
    >$file
    for f in `ls [0-9][0-9][0-9]`;
    do
        echo 'perl ~/my_program3/src/assembly/czl_exonerate_after_blasn.pl -i '$f' -o '$f'.out -q ../Trinity-GG.fasta -t /home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.masked.fasta -p $SLURM_CPUS_PER_TASK -fl 2000' >> $file
    done
    cd ../..
    echo 'cd '$exdir >> $file1
    echo 'if ! [ -f all.out.gz ]; then cat [0-9][0-9][0-9].*out | gzip -c > all.out.gz; fi' >> $file1
    echo 'if ! [ -f all.out.f.gz ]; then ~/my_program3/src/annot_genome/czl_exonerate_after_blasn_filter.pl -i all.out.gz -o all.out.f -qs ../Trinity-GG.tpm1.fasta.fai -qcov 0.75 -score-frac 0.95 -qiden1 0.98 -qiden2 0.95;  gzip all.out.f;  fi' >> $file1
    echo '~/my_program3/src/annot_genome/czl_exonerate_to_gtf_gff3.pl -i all.out.f.gz -o all.out.f' >> $file1
    echo 'cd ../../' >> $file1
done
cmd=r03.ex2.cat; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -g 8 -b 5 --time=1:00:00 -J $cmd --logdir=$cmd.log -f $cmd.sh

file=r03.ex2.cat.fix.sh
>$file;
for sm in `cat sample`
do
    dir=$sm.trinity;
    exdir=$dir/Trinity-GG.carAur01.ex2
#   echo 'cd '$exdir'; cat all.out.f.gtf | awk '"'"'$3!="similarity" && $3!="intron"'"'"' > all.out.f.no_align.gtf' >> $file
    echo 'cd '$exdir'; cat all.out.f.gff3 | awk '"'"'$3!="similarity"'"'"' > all.out.f.no_align.gff' >> $file
done

file=r03.cuffcompare.sh;
>$file;
echo 'cuffcompare -o all.f.ex2.cuffcompare `ls -v *.trinity/Trinity-GG.carAur01.ex2/all.out.f.no_align.gff`;' >> $file

file=r04.filter_mrna.sh
echo "Produce $file"
>$file
IFS=$'\n'
for sm in `cat sample`
do
    dir=$sm.trinity;
    echo 'if ! [ -f '$dir'/tpm1.mrna.tids ]; then
cd '$dir'
>Trinity-GG.mrna_match.ids.tmp
IFS=$'"'"'\n'"'"'
for db in `cat ../nucl_rna_db ../prot_db`; do
    dbname=`echo $db | cut -f 1`;
    dbfile=`echo $db | cut -f 2`;
    ids=`echo $db | cut -f 4`;
    zcat Trinity-GG.$dbname.blastn.m6.gz | cut -f 1,2 | sort -k2,2 -k1,1 > Trinity-GG.$dbname.blastn.m6.id_pairs
    if [ "x$ids" == "x" ]; then
        cut -f 1 Trinity-GG.$dbname.blastn.m6.id_pairs | uniq >> Trinity-GG.mrna_match.ids.tmp; 
    else
        join -t $'"'"'\t'"'"' -1 2 -2 1 Trinity-GG.$dbname.blastn.m6.id_pairs $ids | cut -f 2 | sort | uniq >> Trinity-GG.mrna_match.ids.tmp; 
    fi
    rm Trinity-GG.$dbname.blastn.m6.id_pairs
done
cat Trinity-GG.mrna_match.ids.tmp trinotate_coding_or_long_orf.tids | sort | uniq > Trinity-GG.mrna.tids
cat tpm1.tid Trinity-GG.mrna.tids  | sort | uniq -d > tpm1.mrna.tids
cat tpm1.tid tpm1.mrna.tids | sort | uniq -u > tpm1.ncrna.tids
rm  Trinity-GG.mrna_match.ids.tmp 
cd ../
fi' >> $file
done
cmd=r04.filter_mrna; if ! [ -d $cmd.log ]; then mkdir $cmd.log; else rm $cmd.log/*; fi; swarm -m blast,TransDecoder -g 8 -b 22 --time=1:00:00 -J $cmd --logdir=$cmd.log -f $cmd.sh

file=r05.ex2.filter_mrna.sh
>$file
for sm in `cat sample`
do
    dir=$sm.trinity;
    exdir=$dir/Trinity-GG.carAur01.ex2
    echo 'cd '$exdir >>$file
    echo '~/my_program3/src/annot_genome/czl_exonerate_after_blasn_gff3_remove_dup.pl -i all.out.f.no_align.gff -o all.out.f.no_align.rmdup.gff;'  >> $file
    echo '~/my_program3/src/annot_genome/czl_exonerate_after_blasn_gff3_select_tid.pl -i all.out.f.no_align.rmdup.gff -o all.out.f.no_align.ncrna.gff --list ../tpm1.ncrna.tids' >> $file
    echo 'cd ../../' >> $file
done

file=r06.cuffcompare.sh
>$file
echo 'cuffcompare -o all.f.ex2.ncrna.cuffcompare *.trinity/Trinity-GG.carAur01.ex2/all.out.f.no_align.ncrna.gff' >> $file
echo 'f=all.f.ex2.ncrna.cuffcompare.combined.gtf;' >> $file
echo 'cat $f | perl -e '"'"'
while(<>) {
    if (/^#/ || /^\s*$/) {next; } chomp; @t=split "\t"; 
    foreach $aa (split /\s*;\s*/,$t[8]) {
        my ($u,$v) = split /\s+/, $aa;
        if ($u eq "gene_id") { $gid=$v; }
        elsif ($u eq "geneID" && !defined $gid) { $gid=$v; }
        if ($u eq "transcript_id") { $tid=$v; }
    }
    if ($gid=~m/"(.*)"/) { $gid=$1; }
    if ($tid=~m/"(.*)"/) { $tid=$1; }
    print "$gid\t$tid\t", $_, "\n";
}'"'"' > "$f.XLOC_TCONS"' >> $file
echo 'cat $f | perl -e '"'"'my ($chr0, $b0,$e0,$gid0,$s0);
while(<>) {
    if (/^#/ || /^\s*$/) {next; } chomp; @t=split "\t"; 
    foreach $aa (split /\s*;\s*/,$t[8]) {
        my ($u,$v) = split /\s+/, $aa;
        if ($u eq "gene_id") { $gid=$v; }
        elsif ($u eq "geneID" && !defined $gid) { $gid=$v; }
    }
    if ($gid=~m/"(.*)"/) { $gid=$1; }
    if (defined $gid0 && $gid0 eq $gid) { if ($b0>$t[3]) {$b0=$t[3];} if ($e0<$t[4]) {$e0=$t[4];} }
    else {
        if (defined $gid0) { print "$chr0\t", $b0-1, "\t$e0\t$gid0\t0\t$s0", "\n"; }
        $chr0=$t[0]; $gid0=$gid; $b0=$t[3]; $e0=$t[4]; $s0=$t[6];
    }
}'"'"' > "$f.XLOC.bed"' >> $file
echo 'bedtools subtract -a "$f.XLOC.bed" -b ../carAur01/carAur01.gene.annot.bed -A -f 0.5 | sort -k4,4 > $f.XLOC.not_gene.bed' >> $file 
echo 'join -t$'"'"'\t'"'"' -1 4 -2 1 "$f.XLOC.not_gene.bed" "$f.XLOC_TCONS" | cut -f 8-16 > $f.not_gene.gtf' >> $file
echo 'gffread -K -M $f.not_gene.gtf > $f.not_gene.gff' >> $file
