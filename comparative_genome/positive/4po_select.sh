rm nnn4 


cat nnn3|sed "/^$/d"|sed 's/.*prob = //g'|sed 's/ = //g'|sed 'N;s/\n/\t/g'|awk 'BEGIN{OFS=FS="\t"}{print $0,$2/2}'|awk 'BEGIN{OFS=FS="\t"}{if($3<0.05) print $1}' >nnn4
