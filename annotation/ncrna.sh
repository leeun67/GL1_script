#!/usr/bin/bash
#SBATCH -J ncrna                   # 作业名为 test
#SBATCH -o nc.out               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 24




genome_size_2fold=
genome=
pre=


source /home/a2431559261/OLDISK/soft/miniconda3/etc/profile.d/conda.sh
conda activate ncrna




cmscan -Z ${genome_size_2fold}  --cut_ga --rfam --nohmmonly \
--tblout ${pre}.tblout --fmt 2 --cpu 24  \
 --clanin ~/OLDISK/dataset/Rfam/rfam2/Rfam.clanin \
~/OLDISK/dataset/Rfam/rfam2/Rfam.cm  ${genome} > ${pre}.cmscan


awk 'BEGIN{OFS="\t";}{if(FNR==1) print "target_name\taccession\tquery_name\tquery_start\tquery_end\tstrand\tscore\tEvalue"; if(FNR>2 && $20!="=" && $0!~/^#/) print $2,$3,$4,$10,$11,$12,$17,$18; }' ${pre}.tblout >${pre}.tblout.final.xls


tRNAscan-SE -o ${pre}.tRNA.out -f ${pre}.tRNA.ss -m ${pre}.tRNA.stats ${genome}


barrnap --kingdom euk  --threads 24 --quiet genome.fasta  --outseq ${pre}.rRNA.fasta  > ${pre}.rRNA.gff3

