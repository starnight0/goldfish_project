datadir=/home/chenz11/data/goldfish/11549472/sergey_canu70x/arrow/miRBase_run1
cwd=`pwd`
cd $datadir
echo 0 > id
echo '##gff-version 3' >carAur01.rfam.gff
head -n 2 rfam_by_2Mbp/1.cmscan.tblout > rfam.cmscan.tblout
for i in `seq 1 684`
do
	tail -n +3 rfam_by_2Mbp/$i.cmscan.tblout | grep -v '^#' >> rfam.cmscan.tblout
	# Contig, source, mdl, start, end, 
	tail -n +3 rfam_by_2Mbp/$i.cmscan.tblout | awk -v OFS=$'\t' -v id=`head -n 1 id` '
	{
		if ($14==1) {
			id++;
			id1=sprintf("%08i",id);
			a="ID=RNA"id1":"$3":"$2":"$8"-"$9":"$12":"$17";Name="$3":"$2":"$8"-"$9":"$12";"
			a=a"idx="$1";uRfamID="$3";uRfamName="$2";"
			if ($5 != "-") { a=a"acc="$5";" }
			if ($6 != "-") { a=a"clanID="$6";" }
			a=a"mdl="$7";mdlFrom="$8";mdlTo="$9";trunc="$13";pass="$14";gc="$15";bias="$16";ev="$18";inc="$19";olp="$20";"
			if ($27 != "-") { a=a"RfamDesc="$27";" }
			if ($21 != "-") {
				a=a"anyidx="$21";afrac1="$22";afrac2="$23";"
				if ($24 ~ /[0-9]/) { a=a"winidx="$24";wfrac1="$25";wfrac2="$26";" }
			}
			if ($12=="-") {
				print $4,"cmscan_rfam","mRNA",$11,$10,$17,$12,".",a
			} else {
				print $4,"cmscan_rfam","mRNA",$10,$11,$17,$12,".",a
			}
		}
	} END {print id > "id"}' >> carAur01.rfam.gff
done
head -n 2 rfam.cmscan.tblout > rfam.cmscan.f.tblout
# opl (column 20) not "=", score >=30, E-value < 10e-6
cat rfam.cmscan.tblout | tail -n +3 | awk '$20!="=" && $17>=30 && $18<10e-6 {print}' >> rfam.cmscan.f.tblout
cat rfam.cmscan.f.tblout | tail -n +3 | awk '{print $2}' | sort | uniq -c | awk '{print $2"\t"$1}'
cd $cwd
