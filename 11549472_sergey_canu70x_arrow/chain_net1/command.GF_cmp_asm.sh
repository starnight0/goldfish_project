sp1=carAur01
spl1=~/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.withM.sm.fa.by_100Mbp/files
sp1_2bit=~/data/goldfish/11549472/sergey_canu70x/arrow/carAur01/carAur01.2bit
sp2=fish1p
spl2=/data/projects/burgess/zelin/goldfish/jim_assembly/fish1p/fish1p_pseudohap.1.fasta.by_100Mbp/files
sp2_2bit=/data/projects/burgess/zelin/goldfish/jim_assembly/fish1p/fish1p_pseudohap.1.2bit
~/my_program3/src/assembly/czl_lastz_chain_net_files_array.pl --ref-list $spl1 --query-list $spl2 --ref-2bit $sp1_2bit --query-2bit $sp2_2bit --out-dir $sp1.vs.$sp2 $opt  --lastz-opt '--gapped --ambiguous=n --step=3 --strand=both --masking=10 --maxwordcount=500 --exact=50 --identity=90..100 --format=axt  --allocate:traceback=1024M' --chain-opt '-faQ -faT -minScore=5000 -linearGap=medium' --recovery


sp1=carAur03
spl1=~/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/carAur03.sm.fa.by_chr/files
sp1_2bit=~/data/goldfish/11549472/sergey_canu70x/arrow/carAur03/carAur03.2bit
sp2=fish1p
spl2=/data/projects/burgess/zelin/goldfish/jim_assembly/fish1p/fish1p_pseudohap.1.fasta.by_100Mbp/files
sp2_2bit=/data/projects/burgess/zelin/goldfish/jim_assembly/fish1p/fish1p_pseudohap.1.2bit
~/my_program3/src/assembly/czl_lastz_chain_net_files_array.pl --ref-list $spl1 --query-list $spl2 --ref-2bit $sp1_2bit --query-2bit $sp2_2bit --out-dir $sp1.vs.$sp2 $opt  --lastz-opt '--gapped --ambiguous=n --step=3 --strand=both --masking=10 --maxwordcount=500 --exact=50 --identity=90..100 --format=axt  --allocate:traceback=1024M' --chain-opt '-faQ -faT -minScore=5000 -linearGap=medium' --recovery
