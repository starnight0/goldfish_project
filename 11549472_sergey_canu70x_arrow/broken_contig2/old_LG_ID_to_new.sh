cwd=`pwd`
dwd=`echo $cwd | sed 's/data_cmd/data/'`
idir=$dwd/rmdup.q5_m2.10K/out.LOD_5.5/orig

cat $idir/*.ord | head -n 1 > $idir/all.newid.order
cat old_LG_ID_to_new.tab | perl -ne '
	chomp;
	@t=split "\t";
	if ($t[0] =~ m/^([0-9]+)\.([0-9]+)/) {
		$id0=$1; $sub=$2;
		open IN, "<'$idir'/LG$id0.break.ord";
		while (<IN>) {
			@t1=split "\t";
			if ( ($t1[14] cmp $sub)==0 ) {
				$t1[14]=$t[1];
				print join("\t", @t1);
			}
		}
		close IN;
	} elsif ($t[0] =~ m/^([0-9]+):([0-9]+),([0-9]+)/) {
		$id1=$2; $id2=$3; $id0=$1;
		open IN, "<'$idir'/LG$id0.${id1}_${id2}.ord";
		while (<IN>) {
			@t1=split "\t";
			if ( ($t1[14] cmp $id0)==0 ) {
				$t1[14]=$t[1];
				print join("\t", @t1);
			}
		}
		close IN;
	} elsif ($t[1] ne "0") {
		open IN, "<'$idir'/LG$t[0].order";
		while (<IN>) {
			@t1=split "\t";
			if ( ($t1[14] cmp $t[0])==0 ) {
				$t1[14]=$t[1];
				print join("\t", @t1);
			}
		}
		close IN;
	}
' >> $idir/all.newid.order
tail -n +2 $idir/all.newid.order | awk -F_ '{print $2"_"$3"\t"$4}' | awk -F$'\t' -v OFS=$'\t' '{print $1,$2,$15,$3}' > $idir/all.newid.map4
cat $idir/all.newid.map4 | sed 's/\(\...\)[0-9]\+$/\1/' | sort -k3,3 -k4,4n -k1,1n -k2,2n | awk -F$'\t' 'BEGIN{tid=""; lg="";cm="";} { if ($3!=lg) {print $0} else if ($4!=cm) {print $0} else if ($1!=tid) {print $0} tid=$1; lg=$3; cm=$4}' > $idir/all.newid.f.map4
