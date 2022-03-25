mkdir -p  muscle
cd muscle
ls ../Single_Copy_Orthologue_Sequences >sin_txt

for i in `cat sin_txt`
do    
muscle -in ../Single_Copy_Orthologue_Sequences/${i} -out ${i}.out &
done

sed -i 's/\.fa$//g' sin_txt
cp ../msa.py .

#python msa.py

######进行cds比对

cd ..
mkdir cds_align
cd cds_align
cp ../muscle/sin_txt .
ln -s /home/a2431559261/OLDISK/dataset/genome/max/cds/* .
ls *.cds >cds.name

for i in `cat sin_txt`
do
        rm t1
        cat ../Single_Copy_Orthologue_Sequences/${i}.fa |grep ">"|sed 's/^>//g'  >t1
        for j in `cat cds.name`
        do
                seqtk subseq $j t1 >>${i}.cdd

        done

done


for i in `cat sin_txt`
do
muscle -in ${i}.cdd -out ${i}.cds.out &
done

rm *.cdd

 
