ref_gff=f153.gff3
query_gff=f153.gff3
ref_prefix=f153
query_prefix=f153_bak
ref_cdna=f153.tr
query_cdna=f153.tr
##gff 转换为bed
python -m jcvi.formats.gff bed --type=mRNA  $ref_gff > ${ref_prefix}.bed
python -m jcvi.formats.gff bed --type=mRNA  $query_gff > ${query_prefix}.bed

## 将bed 进行去重复
python -m jcvi.formats.bed uniq ${ref_prefix}.bed
python -m jcvi.formats.bed uniq ${query_prefix}.bed

##得到了ath.uniq.bed和osa.uniq.bed, 根据bed文件第4列就可以用亍提叏cds序列和蛋白序列
seqkit grep -f <(cut -f 4 ${ref_prefix}.uniq.bed ) $ref_cdna | seqkit seq -i > ${ref_prefix}.cds
seqkit grep -f <(cut -f 4 ${ref_prefix}.uniq.bed ) $ref_cdna | seqkit seq -i > ${ref_prefix}.pep
seqkit grep -f <(cut -f 4 ${query_prefix}.uniq.bed ) $query_cdna | seqkit seq -i > ${query_prefix}.cds
seqkit grep -f <(cut -f 4 ${query_prefix}.uniq.bed ) $query_cdna | seqkit seq -i > ${query_prefix}.pep

mkdir -p cds && cd cds
ln -s ../${ref_prefix}.cds .
ln -s ../${ref_prefix}.uniq.bed ${ref_prefix}.bed
ln -s ../${query_prefix}.cds .
ln -s ../${query_prefix}.uniq.bed ${query_prefix}.bed
python -m jcvi.compara.catalog ortholog --no_strip_names ${ref_prefix} ${query_prefix}


python -m jcvi.compara.synteny screen --minspan=30 --simple ${ref_prefix}.${query_prefix}.anchors ${ref_prefix}.${query_prefix}.anchors.new

