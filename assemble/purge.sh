minimap2 -ax map-pb genome.fa subreads.fasta.gz \
    | samtools view -hF 256 - \
    | samtools sort -@ 8 -m 1G -o aligned.bam -T tmp.ali


purge_haplotigs  readhist  -b aligned.bam  -g genome.fasta  -t 20


purge_haplotigs  contigcov  -i aligned.bam.genecov  -o coverage_stats.csv  -l $low_cutoff  -m $mid_cutoff  -h $high_cutoff


purge_haplotigs purge  -g genome.fasta  -c coverage_stats.csv  -b aligned.bam  -t 4  -a 60


