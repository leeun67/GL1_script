#!/usr/bin/bash
#SBATCH -J di                   # 作业名为 test
#SBATCH -o test.out               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 24

export PATH="$PATH:/home/a2431559261/OLDISK/soft/wtdbg2"

gff=FL_H.gff3
pep=fh.pep
cds=fh.cds


#mkdir blast
#mkdir diamond
#cd diamond
#ln -s ../${pep} .

#swiss 比对
mkdir -p blast
cd blast
ln -s ../${pep} .
ln -s ../${cds} .
ln -s ../${gff} .
diamond blastp --db /OLDISK_R730/annotation_database/swissprot/uniprot_sprot.diamond.dmnd --query $pep --outfmt 0 --sensitive  \
 -k 1 --evalue 1e-5 --id 40 --block-size 20.0  --index-chunks 1 --unal 0 --un swiss_unalign.txt --threads 24 --out swiss.txt

cat  swiss_unalign.txt |grep ">" >title_swiss_unalign
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"~":$0 }' ${pep} |grep --file title_swiss_unalign |tr "~" "\n" >swiss_unalign.fasta

#nr 比对
diamond blastp --db /OLDISK_R730/annotation_database/nr/diamond/nr.diamond --query swiss_unalign.fasta --outfmt 0 --sensitive  -k 1 --evalue 1e-5 --id 30 --block-size 20.0  --index-chunks 1 --unal 0 --un nr_unalign.txt --threads 24 --out nr.txt


rm swiss_unalign.fasta swiss_unalign.txt  title_swiss_unalign

#得到blast 简明结果
echo "Query=" >config
echo ">" >>config
cat swiss.txt|grep --file config >swiss.blast
cat nr.txt|grep --file config >nr.blast

#得到swiss blast 综合结果
cat swiss.blast |sed 'N;s/\n/\t/'|cut -f 2|sed 's/ /\t/1'|cut -f 2|sed 's/ OS/\t/' |cut -f 1 >ss
cat swiss.blast |sed 'N;s/\n/\t/'|cut -f 1  >s1
paste s1 ss>s_fi
cat s_fi|sed 's/Query= //g'|sed 's/ gene=/\t/g'|cut -f 1,2 >swiss_result
rm s_fi s1 ss

cat nr.blast |sed 'N;s/\n/\t/'|cut -f 2|sed 's/ /\t/1'|cut -f 2|sed 's/ \[/\t/' |cut -f 1 >ss
cat nr.blast |sed 'N;s/\n/\t/'|cut -f 1  >s1
paste s1 ss>s_fi
cat s_fi|sed 's/Query= //g'|sed 's/ gene=/\t/g'|cut -f 1,2 >nr_result
rm  s1 s_fi ss

cat nr_unalign.txt |grep ">" |cut -d " " -f 1|sed 's/>//g'|awk '{OFS=FS="\t"}{print $1,"Uncharacterized protein"}'>unchar.txt

cat swiss_result nr_result unchar.txt >blast_result



#得到pep注释

cat ${pep}| sed 's/ /\t/1'|cut -f 1|sed 's/^>/>\t/g' >tmp

awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$2}ARGIND==2{if($1==">") print $1,$2,a[$2]; else {print $0}}' blast_result tmp |sed 's/\t//1' >pep.annote
rm tmp

##得到cds注释
cat ${cds}| sed 's/ /\t/1'|cut -f 1|sed 's/^>/>\t/g' >tmp

awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$2}ARGIND==2{if($1==">") print $1,$2,a[$2]; else {print $0}}' blast_result tmp |sed 's/\t//1' >cds.annote
rm tmp



cat $gff |grep -v "#" |sed '/^$/d' > gff

## gff 功能注释
cat blast_result |sed 's/\./\t/1' |awk 'BEGIN{OFS=FS="\t"}{if($2==1) print $0}'|cut -f 1,3 >t.1

cat blast_result |sed 's/\./\t/1' |awk 'BEGIN{OFS=FS="\t"}{if($2==1) print $0}'|cut -f 1 >tt




paste t.1 tt |sed 's/\t/\tNote=/1'|sed 's/^/ID=/g'|sed 's/\t/\tID=/2'>ttt
cat gff|cut -f 3,9|tr ";" "\t" >tmp.gff
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$3]=$2}ARGIND==2{if($1=="gene") print $1,$2,$3,a[$2]; else {print $0}}' ttt tmp.gff >tmp1

cat blast_result|sed 's/\t/\tNote=/1'|sed 's/^/ID=/g'>mt

awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$2}ARGIND==2{if($1=="mRNA") print $1,$2,$3,$4,a[$2]; else {print $0}}' mt tmp1 |cut -f 2,3,4,5,6 |tr "\t" ";" >tmp2
cut -f 1,2,3,4,5,6,7,8 gff >tmp3
paste tmp3 tmp2 >final.gff

rm tt ttt tmp.gff tmp1 tmp2 tmp3 mt 











#cd ..

#cd blast
#ln -s ../${pep} .
#blastp -db /OLDISK_R730/annotation_database/nr/blast/nr.blast -query $pep  -evalue 1e-5 -outfmt 5 -num_threads 20  -out bla_nr.xml> blast_nr.log &

#blastp -db /OLDISK_R730/annotation_database/swissprot/uniprot_sprot.blast -query $pep  -evalue 1e-5 -outfmt 6 -num_threads 20  -out bla_swiss.xml> blast_swill.log &






