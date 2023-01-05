genome_name=goldfish.arrow.renamed
base=$genome_name.maker.output
## fetch finish and fail jobs and contigs
cat $base/${genome_name}_master_datastore_index.log | grep FINISHED | cut -f 1,2 | sort | uniq > contig2.success
cut -f 1 contig2.success | sort | uniq > contig.success

i=0
>tmp.a
for f in `cat ../by_2Mbp.list`
do
	i=`expr $i + 1`
	cat $f | sed -n '/^>/ s/^>//p' | awk '{print $0"\t"'$i'}' >> tmp.a
done
sort -k1,1 tmp.a >tmp.b; mv tmp.b tmp.a
join -j 1 -t $'\t' -v 1 tmp.a contig.success > contig.failed2
cut -f 2 contig.failed2 | sort -k1,1n | uniq > contig.failed.i
fail=`cut -f 2 contig.failed2 | sort -k1,1n | uniq  | head -n 1`
#for i in `cut -f 2 contig.failed2 | sort -k1,1n | uniq | tail -n +2`
#do
#  echo $i
#  fail="$fail,$i"
#done

exit 0
