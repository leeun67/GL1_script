awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$0}ARGIND==2{print $0,a[$1]}' godb.txt ge2ge.tmp |cut -f1,2,3,5,6,7|sed '1 i GOID\tGene\tCount\tDEFINITION\tONTOLOGY\tTERM' >go2gene.txt


awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$0}ARGIND==2{print $1,a[$2]}' godb.txt tmp3 |sed '1 i Gene\tGOID\tDEFINITION\tONTOLOGY\tTERM' >go_annot.txt

cp ~/OLDISK/dataset/genome/2_cb5/go/go_class.txt .

rm ge2ge.tmp
rm tmp3 tmp4



