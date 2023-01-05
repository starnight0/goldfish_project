mkdir detect_het
~/my_program3/src/annot_genome/czl_goldfish_het_gene.pl --bgp big/carAur01.noM.gene.bgp --ctg-info carAur01.ctgInfo2 --self-m6 ../ohnolog/fish17.blastn/pairs.gene.f2/GF.GF.f2.join.m6.gz -o detect_het/ --m6-with-sp --iden1 95

#cut -f 1 detect_het/het_remove_contig.txt > a
#cut -f 14 detect_het/het_remove_contig.txt | sed 's/;/\n/g' | sed 's/,.*$//' >> a
#cat a | sort | uniq > detect_het/tmp.ctg_id
#~/my_program3/src/utility/czl_fasta_fetch -i carAur01.noM.sm.fa -m 1 --list detect_het/tmp.ctg_id -o detect_het/tmp.ctg.
~/my_program3/src/annot_genome/czl_goldfish_het_gene.pl --bgp big/carAur01.noM.gene.bgp --ctg-info carAur01.ctgInfo2 --self-m6 ../ohnolog/fish17.blastn/pairs.gene.f2/GF.GF.f2.join.m6.gz -o detect_het/ --m6-with-sp --iden1 95 --ctg-fa carAur01.noM.sm.fa


mkdir detect_het/run2
~/my_program3/src/annot_genome/czl_goldfish_het_gene.pl --bgp big/carAur01.noM.gene.bgp --ctg-info carAur01.ctgInfo2 --self-m6 ../ohnolog/fish17.blastn/pairs.gene.f2/GF.GF.f2.join.m6.gz --gene-align-stat detect_het/gene_pairs/gene_pairs.stat -o detect_het/run2/ --m6-with-sp --iden1 97.5 --gene-iden 95 


##################################################
# run4 : 2018-12-06
##################################################
# {{{
mkdir detect_het/run4
~/my_program3/src/assembly/czl_detect_het.pl -i ../chain_net/carAur.vs.carAur.I90/all.minScore50000.psl.gz --ctg-info carAur01.ctgInfo2 -o detect_het/run4/  --total_iden_thres1l 0.96 --total_iden_thres2l 0.8  --total_iden_thres1h 0.96  --total_iden_thres2h 0.8 --coverage-thres1 0.5  --coverage-thres2 0.8 2>detect_het/run4/stderr
cat detect_het/run4/pairs.txt | awk -F$'\t' -v OFS=$'\t' '{print $3,$5,$6,$8":"$10":"$11,0,"+"}' | sort -k1,1 -k2,2n > detect_het/run4/pairs.bed
cat detect_het/run4/pairs.txt | awk -F$'\t' -v OFS=$'\t' '{print $8,$10,$11,$3":"$5":"$6,0,$13}' | sort -k1,1 -k2,2n > detect_het/run4/pairs.pri.bed
cat detect_het/run4/pairs.txt | cut -f 3 | sort | uniq > detect_het/run4/het_ctgs
cat detect_het/run4/pairs.txt | cut -f 8 | sort | uniq > detect_het/run4/pri_ctgs
join -t$'\t' -j 1 detect_het/run4/het_ctgs big/carAur01.noM.gene.bgp > detect_het/run4/carAur01.noM.gene.masked.bgp
cat carAur01.noM.sm.fa.fai | cut -f 1 | sort | uniq > carAur01.noM.ctg_ids
cat carAur01.noM.ctg_ids detect_het/run4/het_ctgs | sort | uniq -u > detect_het/run4/nonhet_ctgs 
join -t$'\t' -j 1 detect_het/run4/nonhet_ctgs big/carAur01.noM.gene.bgp > detect_het/run4/carAur01.noM.gene.unmasked.bgp
for t in masked unmasked
do
    cat detect_het/run4/carAur01.noM.gene.$t.bgp | cut -f 4 | sort | uniq > detect_het/run4/carAur01.noM.gene.$t.tids
    cat detect_het/run4/carAur01.noM.gene.$t.bgp | cut -f 18 | sort | uniq > detect_het/run4/carAur01.noM.gene.$t.gids
done
# }}}


