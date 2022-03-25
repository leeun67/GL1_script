#使用说明
#建立三个文件夹
# genes 存放gtf文件
# genome 存放基因组文件
# samples 存放测序数据文件


#软件配置
#conda create -n rnaseq  gffcompare  hisat2 stringtie samtools fsatqc fastp -c bioconda -c conda-forge
#

refgenome=fl.fasta
refgtf=fl.gtf
sample_list=sample_list.txt


#1 原始数据质控
cd samples

#mkdir -p qc
#mkdir -p json_html

#for sample_name in `cat sample_list.txt`
#do
#fastp -i ${sample_name}_1.fastq.gz \
#-I ${sample_name}_2.fastq.gz \
#-o ./qc/${sample_name}_first.fastq.gz \
#-O ./qc/${sample_name}_second.fastq.gz \
#-t 1 \
#-c \
#-p \
#-5 \
#-3  \
#-j ./json_html/${sample_name}.json \
#-h ./json_html/${sample_name}.html
#done

cd ..


#2  hisat2建立索引
mkdir -p index 
cd index
#extract_splice_sites.py ../genes/${refgtf} > ../genes/${refgenome}.ss
#extract_exons.py ../genes/${refgtf} > ../genes/${refgenome}.exon

#hisat2-build   --ss ../genes/${refgenome}.ss --exon ../genes/${refgenome}.exon  ../genome/${refgenome} ${refgenome}
cd ..

#3比对

mkdir -p align
cd align

for sample_name in `cat ../samples/sample_list.txt`
do
hisat2 \
-p 24 \
--dta \
-x ../index/${refgenome} \
-1 ../samples/qc/${sample_name}_first.fastq.gz \
-2 ../samples/qc/${sample_name}_second.fastq.gz \
| samtools sort -@ 24 -o ${sample_name}.bam -  2>${sample_name}.uniq_mapping_rate.txt

samtools index  ${sample_name}.bam

done

cd ..


#4 转录本拼接
mkdir -p expre_count
cd expre_count

for sample_name in `cat ../samples/sample_list.txt`
do
stringtie \
-p 24 \
-G ../genes/${refgtf} \
-o ${sample_name}.gtf \
-l ${sample_name} \
../align/${sample_name}.bam
done



#5 转录本整合与鉴定新转录本
ls -1 *.gtf>mergelist.txt

stringtie --merge \
-p 24 \
-G ../genes/${refgtf} \
-o stringtie_merged.gtf \
mergelist.txt

#转录本比较
gffcompare \
-r ../genes/${refgtf} \
-G \
-o ./gffcompare_results \
stringtie_merged.gtf


#7 重新定量
mkdir -p re_expre_count
cd re_expre_count

for sample_name in `cat ../../samples/sample_list.txt`
do
mkdir ${sample_name}
cd ${sample_name}
stringtie \
-e -B -p 24 \
-G ../../../genes/${refgtf} \
-o ${sample_name}_transcript.gtf  \
../../../align/${sample_name}.bam
cd ..
done

cd ..

#8获取转录本和基因的表达矩阵
python3 /home/a2431559261/OLDISK/script/rna_seq/prepDE.py3 -i re_expre_count/ -g gene_count_matrix.csv -t transcript_count_matrix.csv

cd ..


mkdir deseq2
cd deseq2

cp ../expre_count/gene_count_matrix.csv .
cp ../expre_count/transcript_count_matrix.csv .







