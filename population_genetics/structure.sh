

vcftools --vcf ../fl.pass.snps.genotype.vcf --plink --out final

plink --file final --make-bed --out final


for i in {2..9}
do
    admixture --cv final.bed $i | tee log${i}.out &
done

