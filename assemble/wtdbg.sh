export PATH="$PATH:/home/a2431559261/OLDISK/soft/wtdbg2/"


#进行基因组装
wtdbg2 -t 32 -i test.subreads.bam.fasta -fo fl -L 10000 -g 500m -x sq

#得到一致性序列
wtpoa-cns -t 32 -i fl.ctg.lay.gz -fo fl.raw.fa

#利用三代reads的比对结果对基因组序列进行打磨修正
minimap2 -t 32 -x map-pb -a fl.raw.fa test.subreads.bam.fasta | samtools sort -@ 12 >fl.raw.sort.bam 

samtools view -@ 10 fl.raw.sort.bam | wtpoa-cns -t 32 -d fl.raw.fa -i - -fo fl.cns.fa

