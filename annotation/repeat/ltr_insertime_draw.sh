genome=fl.fasta
title=GL1
cat ${genome}.pass.list |grep -v '#'|cut -f 12 >inser_time


python3 ~/OLDISK/script/repeat/time_draw.py ${title}

