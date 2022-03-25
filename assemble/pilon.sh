bwa index -p index/draft draft.fa
bwa mem -t 20 index/draft read_1.fq.gz read_2.fq.gz | samtools sort -@ 10 -O bam -o align.bam
samtools index -@ 10 align.bam


sambamba markdup -t 10 align.bam align_markdup.bam


samtools view -@ 10 -q 30 align_markdup.bam > align_filter.bam
samtools index -@ 10 align_filter.bam


java -Xmx1000G -jar pilon-1.23.jar --genome draft.fa --frags align_filer.bam \
    --fix snps,indels \
    --output pilon_polished --vcf &> pilon.log
