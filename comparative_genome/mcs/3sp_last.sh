ref_gff=fl.gff3
query1_gff=cb5.gff3
query2_gff=f153.gff3
ref_prefix=fl
query1_prefix=cb5
query2_prefix=f153
ref_cdna=fl.tr
query1_cdna=cb5.tr
query2_cdna=f153.tr

##gff 转换为bed
python -m jcvi.formats.gff bed --type=mRNA  $ref_gff > ${ref_prefix}.bed
python -m jcvi.formats.gff bed --type=mRNA  $query1_gff > ${query1_prefix}.bed
python -m jcvi.formats.gff bed --type=mRNA  $query2_gff > ${query2_prefix}.bed

## 将bed 进行去重复
python -m jcvi.formats.bed uniq ${ref_prefix}.bed
python -m jcvi.formats.bed uniq ${query1_prefix}.bed
python -m jcvi.formats.bed uniq ${query2_prefix}.bed

##得到了ath.uniq.bed和osa.uniq.bed, 根据bed文件第4列就可以用亍提叏cds序列和蛋白序列
seqkit grep -f <(cut -f 4 ${ref_prefix}.uniq.bed ) $ref_cdna | seqkit seq -i > ${ref_prefix}.cds
seqkit grep -f <(cut -f 4 ${ref_prefix}.uniq.bed ) $ref_cdna | seqkit seq -i > ${ref_prefix}.pep
seqkit grep -f <(cut -f 4 ${query1_prefix}.uniq.bed ) $query1_cdna | seqkit seq -i > ${query1_prefix}.cds
seqkit grep -f <(cut -f 4 ${query1_prefix}.uniq.bed ) $query1_cdna | seqkit seq -i > ${query1_prefix}.pep
seqkit grep -f <(cut -f 4 ${query2_prefix}.uniq.bed ) $query2_cdna | seqkit seq -i > ${query2_prefix}.cds
seqkit grep -f <(cut -f 4 ${query2_prefix}.uniq.bed ) $query2_cdna | seqkit seq -i > ${query2_prefix}.pep

mkdir -p cds && cd cds
ln -s ../${ref_prefix}.cds .
ln -s ../${ref_prefix}.uniq.bed ${ref_prefix}.bed
ln -s ../${query1_prefix}.cds .
ln -s ../${query1_prefix}.uniq.bed ${query1_prefix}.bed
ln -s ../${query2_prefix}.cds .
ln -s ../${query2_prefix}.uniq.bed ${query2_prefix}.bed

###find ortholog
python -m jcvi.compara.catalog ortholog --no_strip_names ${query1_prefix} ${ref_prefix}
python -m jcvi.compara.catalog ortholog --no_strip_names ${ref_prefix} ${query2_prefix}


##build sinple
python -m jcvi.compara.synteny screen --minspan=30 --simple ${query1_prefix}.${ref_prefix}.anchors ${query1_prefix}.${ref_prefix}.anchors.new

python -m jcvi.compara.synteny screen --minspan=30 --simple ${ref_prefix}.${query2_prefix}.anchors ${ref_prefix}.${query2_prefix}.anchors.new






