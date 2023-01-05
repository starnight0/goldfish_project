datadir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/miRBase_run1
cwd=`pwd`

genome1fa=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/goldfish.arrow.renamed.masked.fasta

cd $datadir

head -n 2 rfam.cmscan.tblout > rfam.cmscan.f.tblout
# opl (column 20) not "=", score >=30, E-value < 10e-6
cat rfam.cmscan.tblout | tail -n +3 | awk '$20!="=" && $17>=30 && $18<10e-6 {print}' >> rfam.cmscan.f.tblout
cat rfam.cmscan.f.tblout | tail -n +3 | awk '{print $2}' | sort | uniq -c | awk '{print $2"\t"$1}'
~/my_program3/src/annot_genome/czl_cmscan_to_gff.pl -i rfam.cmscan.f.tblout -o carAur01.rfam.cmscan.f.gff
~/my_program3/src/annot_genome/czl_cmscan_to_gff.pl -i rfam.cmscan.f.tblout -o carAur01.rfam.cmscan.f.bgp -of bgp
bedtools sort -i carAur01.rfam.cmscan.f.bgp  > a; mv a carAur01.rfam.cmscan.f.bgp 
bedToBigBed -extraIndex=name,name2,geneName,geneName2 -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as carAur01.rfam.cmscan.f.bgp $genome1fa.fai carAur01.rfam.cmscan.f.bgp.bb


fn=carAur01.rfam
gff3ToGenePred -geneNameAttr=Name -attrsOut=$fn.gff.attrs $fn.gff stdout | \
		genePredToBigGenePred stdin stdout | sort -k1,1 -k2,2n > $fn.bgp
bedToBigBed -extraIndex=name,name2 -type=bed12+8 -tab -as=/data/genome/jksrc_v352/kent/src/hg/lib/bigGenePred.as $fn.bgp $genome1fa.fai $fn.bb 



cd $cwd
