win1=cb5
win2=f153.fasta_merge_win
pre=f153
gff3=../f153_cb5/data/f153.gff3
bedtools intersect -a ${win1} -b ${win2}|awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$3-$2+1}'  >${pre}.common_specwin


wc ${pre}.common_specwin >sta

awk 'BEGIN{total=0}{total+=$4}END{print total}' ${pre}.common_specwin >>sta


python ~/OLDISK/script/sv/pav_gene.py ${pre}.common_specwin  ${gff3}
