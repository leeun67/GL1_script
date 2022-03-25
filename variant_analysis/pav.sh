#!/usr/bin/bash
#SBATCH -J pav_gene                   # 作业名为 test
#SBATCH -o pav.out               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 24


## 前期准备
## 1 准备data文件夹，里面放query基因组和 ref 基因组 的 fasta文件和gff3 注释文件,并构建bwa索引文件
## 2 检查序列文件和gff三文件， 染色体名字不要含有 ~，然后染色体名字中间不要有空格
## 3 gff文件中，基因名字除了末尾，中间不要含有' . '




t=24
query=cxn.fasta
ref=fipak.fasta

query_index_dir=/home/a2431559261/OLDISK/caixin_sv/cnx_pak/data
ref_index_dir=/home/a2431559261/OLDISK/caixin_sv/cnx_pak/data


query_gff3=cxn.gff3
ref_gff3=fipak.gff3





source /home/a2431559261/OLDISK/soft/miniconda3/etc/profile.d/conda.sh
conda activate sv

mkdir -p query
mkdir -p ref

##query
cd query
##1 窗口切割以及比对


awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"%":$0 }' ../data/${query} |tr '%' '\n' > ../data/2${query}
python  ~/OLDISK/script/sv/win.py ../data/2${query} >step_${query}
rm ../data/2${query}


srun -c 24 bwa mem -w 500 -M ${ref_index_dir}/${ref} step_${query} -t $t  >${query}.sam


cat ${query}.sam|grep -v '@' |grep -v '#'|cut -f 1,2,6 >t1

cut -f 1,2 t1 >t2
cut -f 3 t1 |sed 's/M/M\t/g'|sed 's/S/S\t/g'|sed 's/D/D\t/g'|sed 's/P/P\t/g'|sed 's/H/H\t/g'|sed 's/I/I\t/g'|sed 's/N/N\t/g'|sed 's/X/X\t/g' |grep -v '\[' >t3

paste t2 t3 >tmp
rm t1 t2 t3



##2 提取所有窗口的序列名
cat step_${query} |grep '>'|sed 's/>//g' >all_seq


##3 运行PY 文件 得到特异性区域
python ~/OLDISK/script/sv/tj.py  tmp  query_conserved_seq  all_seq
#结果 specific_seq.txt merge_win sta query_conserved_seq


##4 得到PAV基因
python ~/OLDISK/script/sv/pav_gene.py merge_win ../data/${query_gff3}

cat pav.gene|cut -f1|sort|uniq >${query}_pav.gene

wc ${query}_pav.gene >>sta

mv merge_win ${query}_merge_win
rm pav.gene 
mv sta query_sta
rm tmp
rm all_seq
rm 2tmp_cds.bed

cd ../ref
######ref


##1 窗口切割以及比对
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"%":$0 }' ../data/${ref} |tr '%' '\n' > ../data/2${ref}


python  ~/OLDISK/script/sv/win.py ../data/2${ref} >step_${ref}
rm ../data/2${ref}

bwa mem -w 500 -M ${query_index_dir}/${query} step_${ref} -t $t  >${ref}.sam

#cat ${ref}.sam|grep -v '@' |grep -v '#'|cut -f 1,2,6|sed 's/M/M\t/g'|sed 's/S/S\t/g'|sed 's/D/D\t/g'|sed 's/P/P\t/g'|sed 's/H/H\t/g'|sed 's/I/I\t/g'|sed 's/N/N\t/g'|sed 's/X/X\t/g'|grep -v '\[' >tmp

cat ${ref}.sam|grep -v '@' |grep -v '#'|cut -f 1,2,6 >t1

cut -f 1,2 t1 >t2
cut -f 3 t1 |sed 's/M/M\t/g'|sed 's/S/S\t/g'|sed 's/D/D\t/g'|sed 's/P/P\t/g'|sed 's/H/H\t/g'|sed 's/I/I\t/g'|sed 's/N/N\t/g'|sed 's/X/X\t/g' |grep -v '\[' >t3

paste t2 t3 >tmp
rm t1 t2 t3




##2 提取所有窗口的序列名
cat step_${ref} |grep '>'|sed 's/>//g' >all_seq


##3 运行PY 文件 得到特异性区域
python ~/OLDISK/script/sv/tj.py  tmp  ref_conserved_seq  all_seq
#结果 specific_seq.txt merge_win sta query_conserved_seq


##4 得到PAV基因
python ~/OLDISK/script/sv/pav_gene.py merge_win ../data/${ref_gff3}

cat pav.gene|cut -f1|sort|uniq >${ref}_pav.gene

wc ${ref}_pav.gene >>sta


mv merge_win ${ref}_merge_win
rm pav.gene  
mv sta ref_sta
rm tmp
rm all_seq
rm 2tmp_cds.bed






