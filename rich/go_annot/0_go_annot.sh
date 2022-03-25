
tsv_file=cb5.pep.tsv

cp ~/OLDISK/dataset/genome/2_cb5/go/godb.txt .
cut -f 1,14 ${tsv_file}|awk '{OFS=FS="\t"}{if($2!="") print $0}'|tr '|' '\t' >tmp

python /home/a2431559261/OLDISK/script/rich/go_annot/1_go_annot.py



