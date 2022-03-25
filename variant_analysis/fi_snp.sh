cat snp.cb5.fasta |cut -f 1,2,3,11|tail -n +5|grep -v '\.'|awk 'BEGIN{OFS=FS="\t"}{print $4,$1}'  >snp


cat snp.cb5.fasta |cut -f 1,2,3,11|tail -n +5|grep  '\.'|awk 'BEGIN{OFS=FS="\t"}{print $4,$1}' >indel


python ~/OLDISK/script/sv/snp_tj.py snp fl.pass.snps.genotype.vcf indel fl.pass.indel.genotype.vcf
