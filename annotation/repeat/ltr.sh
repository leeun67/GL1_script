#!/bin/sh
#SBATCH -J ltr                   # 作业名为 test
#SBATCH -o test.out               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 24



genome=cb5.fasta
pre=cb5



#workdir= `pwd`

#cd /mnt/memorydisk/
#ln -s ${workdir}/${genome} .

#LTRharvest
#~/OLDISK/soft/genometools-1.6.0/bin/gt suffixerator \
#  -db ${genome} \
#  -indexname ${pre} \
#  -tis -suf -lcp -des -ssp -sds -dna
#~/OLDISK/soft/genometools-1.6.0/bin/gt ltrharvest \
#  -index ${pre} \
#  -similar 90 -vic 10 -seed 20 -seqids yes \
#  -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 \
#  -motif TGCA -motifmis 1  > ${pre}.harvest.scn 
# LTR_FINDER
#ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 ${genome} > ${pre}.finder.scn 

LTR_retriever -genome $genome -inharvest ${pre}.harvest.scn -infinder ${pre}.finder.scn -threads 24
#rm ${genome}

#rsync  -auP * ${work_dir}
#rm * -rf

 

#LTR_retriever -genome $genome -inharvest ${pre}.harvest.scn -infinder ${pre}.finder.scn -threads 10

