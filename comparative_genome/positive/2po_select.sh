rm nnn2
rm nnn
for i in  `cat a`
do 
	echo $i >>nnn
	cat ${i}.null.out |grep 'lnL'|sed 's/^.*-//g'|sed 's/  *.*//g' >>nnn
	cat  ${i}.alt.out |grep 'lnL'|sed 's/^.*-//g'|sed 's/  *.*//g' >>nnn
done

cat nnn |sed 'N;N;s/\n/\t/g'|awk 'BEGIN{OFS=FS="\t"}{print $0,$3-$2+$3-$2}'|cut -f 1,4 >nnn2


