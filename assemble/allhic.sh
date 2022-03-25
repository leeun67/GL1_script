#$ git clone https://github.com/tangerzhang/ALLHiC
#$ cd ALLHiC
#$ chmod +x bin/*
#$ chmod +x scripts/*  
#$ export PATH=/your/path/to/ALLHiC/scripts/:/your/path/to/ALLHiC/bin/:$PATH
#conda deactivate
#conda deactivate
#conda activate falcon


export PATH="$PATH:/home/a2431559261/OLDISK/soft/ALLHiC/scripts"
export PATH="$PATH:/home/a2431559261/OLDISK/soft/ALLHiC/bin"
export PATH="$PATH:/home/a2431559261/OLDISK/soft/MCScanX/downstream_analyses/"
export PATH="$PATH:/home/a2431559261/OLDISK/soft/khaper/Bin"

export PATH="$PATH:/home/a2431559261/OLDISK/soft/MCScanX/"
#1,建立索引
samtools faidx draft.asm.fasta 
bwa index -a bwtsw draft.asm.fasta  


#2，序列回帖
bwa aln -t 24 draft.asm.fasta P-G_1.fq.gz > reads_R1.sai  
bwa aln -t 24 draft.asm.fasta P-G_2.fq.gz > reads_R2.sai  
bwa sampe draft.asm.fasta reads_R1.sai reads_R2.sai P-G_1.fq.gz P-G_2.fq.gz > sample.bwa_aln.sam  


perl /home/a2431559261/OLDISK/soft/ALLHiC/scripts/proces_sam_2.pl sample.bwa_aln.sam draft.asm.fasta MBOI
# 如果已有BAM文件
# PreprocessSAMs.pl sample.bwa_aln.bam draft.asm.fasta MBOI
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam
samtools view -bt draft.asm.fasta.fai sample.clean.sam > sample.clean.bam

ALLHiC_partition -b sample.clean.bam -r draft.asm.fasta -e GATC -k 25



ALLHiC_rescue -b sample.clean.bam -r draft.asm.fasta \
    -c sample.clean.clusters.txt \
    -i sample.clean.counts_GATC.txt



allhic extract sample.clean.bam draft.asm.fasta --RE GATC

# 优化
for i in group*.txt; do    allhic optimize $i sample.clean.clm;done


ALLHiC_build draft.asm.fasta  


samtools faidx groups.asm.fasta
cut -f 1,2 groups.asm.fasta.fai  > chrn.list
ALLHiC_plot sample.clean.bam groups.agp chrn.list 500k pdf


source /home/a2431559261/OLDISK/soft/miniconda3/etc/profile.d/conda.sh
ln -s /home/a2431559261/OLDISK/genome_assemble/mcs/cb.dna.fa .

conda activate mcscan
ln -s /home/a2431559261/OLDISK/genome_assemble/mcs/f153_cb5/S1_R1.fastq.gz .
ln -s /home/a2431559261/OLDISK/genome_assemble/mcs/f153_cb5/S1_R2.fastq.gz .

cp /home/a2431559261/OLDISK/script/transcript_two_genome_MCscan.py  .
cp /home/a2431559261/OLDISK/script/mykit.py  .

mv groups.asm.fasta k25.dna.fa
python3 transcript_two_genome_MCscan.py cb k25 S1 d 24

