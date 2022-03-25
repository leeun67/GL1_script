vcftools --vcf final.snp.vcf --plink --out final

plink --file final --pca 3 --out final

cat head final.eigenvec > final.pop.eigenvec

Rscript pca.R final.pop.eigenvec
Rscript pca.R final.group.eigenvec
