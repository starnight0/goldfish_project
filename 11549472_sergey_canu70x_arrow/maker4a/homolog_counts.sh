file=stat.txt
>$file
cd p.blast/sp5.pairs.2
if ! [ -f tmp.ENSDAR.carAur.pair ]; then 
zcat ENSDAR.carAur.f4.join.m6.gz | awk '$1~/ENSDAR/ {print $1"\t"$2}' >  tmp.ENSDAR.carAur.pair;
zcat ENSDAR.carAur.f4.join.m6.gz | awk '$2~/ENSDAR/ {print $2"\t"$1}' >> tmp.ENSDAR.carAur.pair;
fi
if ! [ -f tmp.ENSDAR.CYPCAR.pair ]; then
zcat CYPCAR.ENSDAR.f4.join.m6.gz | awk '$1~/ENSDAR/ {print $1"\t"$2}' >  tmp.ENSDAR.CYPCAR.pair;
zcat CYPCAR.ENSDAR.f4.join.m6.gz | awk '$2~/ENSDAR/ {print $2"\t"$1}' >> tmp.ENSDAR.CYPCAR.pair;
fi
if ! [ -f tmp.carAur.CYPCAR.pair ]; then
zcat CYPCAR.carAur.f4.join.m6.gz | awk '$1~/carAur/ {print $1"\t"$2}' >  tmp.carAur.CYPCAR.pair;
zcat CYPCAR.carAur.f4.join.m6.gz | awk '$2~/carAur/ {print $2"\t"$1}' >> tmp.carAur.CYPCAR.pair;
fi
# N GF genes in ZF and CC
cut -f 2 tmp.ENSDAR.carAur.pair | sort | uniq > ZF_GF.GF
cut -f 1 tmp.ENSDAR.carAur.pair | sort | uniq > ZF_GF.ZF
cut -f 1 tmp.carAur.CYPCAR.pair | sort | uniq > GF_CC.GF
cut -f 2 tmp.carAur.CYPCAR.pair | sort | uniq > GF_CC.CC
cut -f 2 tmp.ENSDAR.CYPCAR.pair | sort | uniq > ZF_CC.CC
cut -f 1 tmp.ENSDAR.CYPCAR.pair | sort | uniq > ZF_CC.ZF
cat ZF_GF.GF tmp.GF_CC.GF | sort | uniq -d > ZF_GF_CC.GF
cat ZF_GF.ZF tmp.ZF_CC.ZF | sort | uniq -d > ZF_GF_CC.ZF
cat ZF_CC.CC tmp.GF_CC.CC | sort | uniq -d > ZF_GF_CC.CC
cat ZF_GF.GF ZF_GF_CC.GF | sort | uniq -u > ZF_GF_only.GF 
cat GF_CC.GF ZF_GF_CC.GF | sort | uniq -u > GF_CC_only.GF 
cat ZF_CC.CC ZF_GF_CC.CC | sort | uniq -u > ZF_CC_only.CC
cat GF_CC.CC ZF_GF_CC.CC | sort | uniq -u > GF_CC_only.CC
cat ZF_GF.ZF ZF_GF_CC.ZF | sort | uniq -u > ZF_GF_only.ZF 
cat ZF_CC.ZF ZF_GF_CC.ZF | sort | uniq -u > ZF_CC_only.ZF 

n1=`wc -l ZF_GF_CC.GF | cut -d" " -f 1`
n2=`wc -l ZF_GF_CC.CC | cut -d" " -f 1`
n3=`wc -l ZF_GF_CC.ZF | cut -d" " -f 1`
echo "ZF_GF_CC shared genes" >> $file
echo "GF"$'\t'"CC"$'\t'"ZF" >> $file
echo $n1$'\t'$n2$'\t'$n3 >> $file
echo "GF genes shared only in two species" >> $file
n1=`wc -l ZF_GF_only.GF | cut -d" " -f 1`
n2=`wc -l GF_CC_only.GF | cut -d" " -f 1`
echo "GF_CC"$'\t'"GF_ZF" >> $file
echo $n1$'\t'$n2 >> $file
n1=`wc -l ZF_CC_only.CC | cut -d" " -f 1`
n2=`wc -l GF_CC_only.CC | cut -d" " -f 1`
echo "CC_GF"$'\t'"CC_ZF" >> $file
echo $n1$'\t'$n2 >> $file
n1=`wc -l ZF_GF_only.ZF | cut -d" " -f 1`
n2=`wc -l ZF_CC_only.ZF | cut -d" " -f 1`
echo "ZF_GF"$'\t'"ZF_CC" >> $file
echo $n1$'\t'$n2 >> $file
cd ../../
