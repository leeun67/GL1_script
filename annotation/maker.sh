#!/usr/bin/bash
#SBATCH -J man                   # 作业名为 test
#SBATCH -o test.out               # 屏幕上的输出文件重定向到 test.out
#SBATCH -p comp                   # 作业提交的分区
#SBATCH --qos=normal               # 作业使用的 QoS 为 debug ,normal。debug优先级高，可使用资源少
#SBATCH -N 1                      # 作业申请几个节点
#SBATCH -c 48


cp /home/a2431559261/OLDISK/soft/gmes_linux_64/.gm_key .

genome=man.dna
pre=man
est=Trinity.fasta
protein=a.pep

source /home/a2431559261/OLDISK/soft/miniconda3/etc/profile.d/conda.sh
conda activate maker

mkdir -p repeat
cd repeat
/home/a2431559261/OLDISK/soft/RepeatModeler-2.0.1/BuildDatabase -name ${pre}db ../data/$genome
/home/a2431559261/OLDISK/soft/RepeatModeler-2.0.1/RepeatModeler -database ${pre}db  -pa 24  -LTRStruct  -ninja_dir /home/a2431559261/OLDISK/soft/NINJA-0.95-cluster_only/NINJA >run.out


mkdir mask1
cd mask1
ln -s ../../data/${genome} .
RepeatMasker -engine ncbi -parallel 24  -nolow -no_is -norna  -lib ../${pre}db-families.fa $genome 

cd ..
mkdir mask2
cd mask2
ln -s ../mask1/${genome}.masked .
RepeatMasker -engine ncbi -parallel 24  -noint -small -no_is -norna  -lib ../${pre}db-families.fa $genome.masked
cd ..
cd ..
cd data
ln -s ../repeat/mask2/${genome}.masked.masked ${genome}.softmask
cd ..

est2genome=1
protein2genome=1
cpus=2
mpiexec=20
always_complete=0
model_org=

#####maker 第一轮
maker --CTL
mv maker_opts.ctl rnd1_maker_opts.ctl
sed -i "s/^model_org=all /model_org=${model_org} /g" rnd1_maker_opts.ctl
sed -i "s/^genome= /genome=data\/${genome}.softmask /g" rnd1_maker_opts.ctl
sed -i "s/^est= /est=data\/${est} /g" rnd1_maker_opts.ctl
sed -i "s/^protein= /protein=data\/${protein} /g" rnd1_maker_opts.ctl
######sed -i "s/^rmlib= /rmlib=data\/${rmlib} /g" ${rnd_num}_maker_opts.ctl

sed -i "s/^est2genome=0 /est2genome=${est2genome} /g" rnd1_maker_opts.ctl
sed -i "s/^protein2genome=0 /protein2genome=${protein2genome} /g" rnd1_maker_opts.ctl
sed -i "s/^cpus=1 /cpus=${cpus} /g" rnd1_maker_opts.ctl
sed -i "s/^always_complete=0 /always_complete=${always_complete} /g" rnd1_maker_opts.ctl
######sed -i "s/^model_org=all /model_org=${model_org} /g" rnd1_maker_opts.ctl
#sed -i "s/RepeatMasker=\/home\/a2431559261\/OLDISK\/soft\/miniconda3\/envs\/maker\/bin\/RepeatMasker/RepeatMasker=\/home\/a2431559261\/OLDISK\/soft\/RepeatMasker\/RepeatMasker/g" maker_exe.ctl

mpiexec -n ${mpiexec} maker -base round1 rnd1_maker_opts.ctl maker_bopts.ctl maker_exe.ctl


cd round1.maker.output
gff3_merge -s -d round1_master_datastore_index.log > rnd1.all.maker.gff
fasta_merge -d round1_master_datastore_index.log
###### GFF w/o the sequences
gff3_merge -n -s -d round1_master_datastore_index.log > rnd1.all.maker.noseq.gff
gff3_merge -g -s -n -d round1_master_datastore_index.log > rnd1.gene_model.maker.noseq.gff
cat rnd1.gene_model.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }' >rnd1_genesta.txt
cat round1.all.maker.proteins.fasta |grep ">"|cut -d " " -f 3|tr ":" "\t"|awk '{print $0 }' | awk '{if($2<0.5) sum+=1}END{ print NR, sum / NR }' >>rnd1_genesta.txt

conda activate busco
run_BUSCO.py -i round1.all.maker.proteins.fasta  -o busco_rnd1 -l /home/a2431559261/OLDISK/dataset/embryophyta_odb10/ -m prot -c 24 -sp arabidopsis
cd ..

export PATH="$PATH:/home/a2431559261/OLDISK/soft/gmes_linux_64/ProtHint/bin"
export PATH="$PATH:/home/a2431559261/OLDISK/soft/gmes_linux_64/"


conda activate genemark

#module load genemark/genemark-4.59
mkdir -p predict
cd predict
#mkdir gm1
#cd gm1
#ln -s ../../round1.maker.output/round1.all.maker.proteins.fasta .
#prothint.py --threads 26 --workdir ./g1 ../../data/${genome}.softmask round1.all.maker.proteins.fasta 
#gmes_petap.pl --EP ./prothint.gff --evidence ./g1/evidence.gff  --cores 26 --max_intron 10000 --seq ../../data/${genome}.softmask  --verbose
#gtf2gff3.pl genemark.gtf > genemark.gff3
#cd ../


conda activate maker

###2.2 SNAP
mkdir -p snap1
cd snap1

###使用置信基因模型
maker2zff -x 0.25 -l 50 -d ../../round1.maker.output/round1_master_datastore_index.log
# -x max AED (0.5 default)
#-l min length of protein sequence
mv genome.dna rnd1.zff.length50_aed0.25.dna
mv genome.ann rnd1.zff.length50_aed0.25.ann

#2  gather some stats and validate 做一些统计和验证
fathom rnd1.zff.length50_aed0.25.ann rnd1.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom rnd1.zff.length50_aed0.25.ann rnd1.zff.length50_aed0.25.dna -validate > validate.log 2>&1

#3 收集训练序列
fathom rnd1.zff.length50_aed0.25.ann rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# -category 将收集的模型基因进行分类，有uni wrn alt err olp,一般只有uni 可用于下一步分析
# -export 导出四个文件  -plus 将导出来的序列转换为正链  具体不懂

#4 create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..

# assembly the HMM
hmm-assembler.pl rnd1.zff.length50_aed0.25 params > rnd1.zff.length50_aed0.25.hmm
cd ../..


###maker 第二轮
#1，回收第一轮的证据
# transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' round1.maker.output/rnd1.all.maker.noseq.gff > rnd1.all.maker.est2genome.gff
# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' round1.maker.output/rnd1.all.maker.noseq.gff > rnd1.all.maker.protein2genome.gff
# repeat alignments
#awk '{ if ($2 ~ "repeat") print $0 }' round1.maker.output/rnd1.all.maker.noseq.gff > rnd1.all.maker.repeats.gff


est=
protein=
rmlib=
est2genome=0
protein2genome=0
est_gff=rnd1.all.maker.est2genome.gff
protein_gff=rnd1.all.maker.protein2genome.gff
#rm_gff=rnd1.all.maker.repeats.gff
trna=0
snaphmm=rnd1.zff.length50_aed0.25.hmm
cpus=2
mpiexec=24
always_complete=1
model_org=
min_protein=30
single_exon=1 ########turn it on for fungi genome annotation
single_length=150

maker -OPT
mv maker_opts.ctl rnd2_maker_opts.ctl


sed -i "s/^min_protein=0 /min_protein=30 /g" rnd2_maker_opts.ctl
#sed -i "s/^single_exon=0 /single_exon=1 /g" rnd2_maker_opts.ctl
#sed -i "s/^single_length=250 /single_length=150 /g" rnd2_maker_opts.ctl

sed -i "s/^est_gff= /est_gff=${est_gff} /g" rnd2_maker_opts.ctl
sed -i "s/^protein_gff= /protein_gff=${protein_gff} /g" rnd2_maker_opts.ctl
#####sed -i "s/^rm_gff= /rm_gff=${rm_gff} /g" rnd2_maker_opts.ctl
#####sed -i "s/^trna=0 /trna=${trna} /g" rnd2_maker_opts.ctl
sed -i "s/^snaphmm= /snaphmm=predict\/snap1\/${snaphmm} /g" rnd2_maker_opts.ctl
#####sed -i "s/^gmhmm= /gmhmm=predict\/gm1\/g1\/GeneMark_ES\/gmhmm.mod /g" rnd2_maker_opts.ctl
#######sed -i "s/^augustus_species= /augustus_species=${augustus_species} /g" ${rnd_num}_maker_opts.ctl

sed -i "s/^genome= /genome=data\/${genome} /g" rnd2_maker_opts.ctl
sed -i "s/^est= /est=${est} /g" rnd2_maker_opts.ctl
sed -i "s/^protein= /protein=${protein} /g" rnd2_maker_opts.ctl
######sed -i "s/^rmlib= /rmlib=${rmlib} /g" ${rnd_num}_maker_opts.ctl
sed -i "s/^est2genome=0 /est2genome=${est2genome} /g" rnd2_maker_opts.ctl
sed -i "s/^protein2genome=0 /protein2genome=${protein2genome} /g" rnd2_maker_opts.ctl
sed -i "s/^cpus=1 /cpus=${cpus} /g" rnd2_maker_opts.ctl
sed -i "s/^always_complete=0 /always_complete=${always_complete} /g" rnd2_maker_opts.ctl
sed -i "s/^model_org=all /model_org=${model_org} /g" rnd2_maker_opts.ctl
sed -i "s/^max_dna_len=100000 /max_dna_len=200000 /g" rnd2_maker_opts.ctl
#####sed -i "s/^split_hit=10000 /split_hit=10000 /g" rnd2_maker_opts.ctl

mpiexec -n ${mpiexec} maker -base round2 rnd2_maker_opts.ctl maker_bopts.ctl maker_exe.ctl



###第三轮
cd round2.maker.output
gff3_merge -s -d round2_master_datastore_index.log > rnd2.all.maker.gff
fasta_merge -d round2_master_datastore_index.log
gff3_merge -n -s -d round2_master_datastore_index.log > rnd2.all.maker.noseq.gff
gff3_merge -g -s -n -d round2_master_datastore_index.log > rnd2.gene_model.maker.noseq.gff
cat rnd2.gene_model.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }' >rnd2_genesta.txt
cat round2.all.maker.proteins.fasta |grep ">"|cut -d " " -f 3|tr ":" "\t"|awk '{print $0 }' | awk '{if($2<0.5) sum+=1}END{ print NR, sum / NR }' >>rnd2_genesta.txt

conda activate busco
run_BUSCO.py -i round2.all.maker.proteins.fasta  -o busco_rnd2 -l /home/a2431559261/OLDISK/dataset/embryophyta_odb10/ -m prot -c 24 -sp arabidopsis
cd ..



conda activate genemark
cd predict
mkdir gm2
cd gm2
ln -s ../../round2.maker.output/round2.all.maker.proteins.fasta .
prothint.py --threads 26 --workdir ./g2 ../../data/${genome}.softmask round2.all.maker.proteins.fasta
gmes_petap.pl --EP ./prothint.gff --evidence ./g2/evidence.gff  --cores 26 --max_intron 10000 --seq ../../data/${genome}.softmask  --verbose
cd ..

conda activate maker
mkdir -p snap2
cd snap2

maker2zff -x 0.25 -l 50 -d ../../round2.maker.output/round2_master_datastore_index.log
mv genome.dna rnd2.zff.length50_aed0.25.dna
mv genome.ann rnd2.zff.length50_aed0.25.ann

fathom rnd2.zff.length50_aed0.25.ann rnd2.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom rnd2.zff.length50_aed0.25.ann rnd2.zff.length50_aed0.25.dna -validate > validate.log 2>&1
fathom rnd2.zff.length50_aed0.25.ann rnd2.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..

hmm-assembler.pl rnd2.zff.length50_aed0.25 params > rnd2.zff.length50_aed0.25.hmm
cd ../..

est=
protein=
rmlib=
est2genome=0
protein2genome=0
#est_gff=rnd2.all.maker.est2genome.gff
#protein_gff=rnd2.all.maker.protein2genome.gff
#rm_gff=rnd1.all.maker.repeats.gff
trna=1
snaphmm=rnd2.zff.length50_aed0.25.hmm
cpus=2
mpiexec=24
always_complete=1
model_org=
min_protein=30
single_exon=1 ########turn it on for fungi genome annotation
single_length=150


maker -OPT
mv maker_opts.ctl rnd3_maker_opts.ctl

sed -i "s/^min_protein=0 /min_protein=30 /g" rnd3_maker_opts.ctl
#sed -i "s/^single_exon=0 /single_exon=1 /g" rnd3_maker_opts.ctl
#sed -i "s/^single_length=250 /single_length=150 /g" rnd3_maker_opts.ctl

sed -i "s/^est_gff= /est_gff=${est_gff} /g" rnd3_maker_opts.ctl
sed -i "s/^protein_gff= /protein_gff=${protein_gff} /g" rnd3_maker_opts.ctl
#####sed -i "s/^trna=0 /trna=${trna} /g" rnd3_maker_opts.ctl
sed -i "s/^snaphmm= /snaphmm=predict\/snap2\/${snaphmm} /g" rnd3_maker_opts.ctl
sed -i "s/^gmhmm= /gmhmm=predict\/gm2\/g2\/GeneMark_ES\/gmhmm.mod /g" rnd3_maker_opts.ctl
sed -i "s/^genome= /genome=data\/${genome} /g" rnd3_maker_opts.ctl
sed -i "s/^est= /est=${est} /g" rnd3_maker_opts.ctl
sed -i "s/^protein= /protein=${protein} /g" rnd3_maker_opts.ctl
sed -i "s/^est2genome=0 /est2genome=${est2genome} /g" rnd3_maker_opts.ctl
sed -i "s/^protein2genome=0 /protein2genome=${protein2genome} /g" rnd3_maker_opts.ctl
sed -i "s/^cpus=1 /cpus=${cpus} /g" rnd3_maker_opts.ctl

sed -i "s/^always_complete=0 /always_complete=${always_complete} /g" rnd3_maker_opts.ctl

sed -i "s/^model_org=all /model_org=${model_org} /g" rnd3_maker_opts.ctl
sed -i "s/^max_dna_len=100000 /max_dna_len=200000 /g" rnd3_maker_opts.ctl
#####sed -i "s/^split_hit=10000 /split_hit=15000 /g" rnd3_maker_opts.ctl
mpiexec -n ${mpiexec} maker -base round3 rnd3_maker_opts.ctl maker_bopts.ctl maker_exe.ctl




#第四轮
cd round3.maker.output
gff3_merge -s -d round3_master_datastore_index.log > rnd3.all.maker.gff
fasta_merge -d round3_master_datastore_index.log
gff3_merge -n -s -d round3_master_datastore_index.log > rnd3.all.maker.noseq.gff
gff3_merge -g -s -n -d round3_master_datastore_index.log > rnd3.gene_model.maker.noseq.gff
cat rnd3.gene_model.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }' >rnd3_genesta.txt
cat round3.all.maker.proteins.fasta |grep ">"|cut -d " " -f 3|tr ":" "\t"|awk '{print $0 }' | awk '{if($2<0.5) sum+=1}END{ print NR, sum / NR }' >>rnd3_genesta.txt

conda activate busco
run_BUSCO.py -i round3.all.maker.proteins.fasta  -o busco_rnd3 -l /home/a2431559261/OLDISK/dataset/embryophyta_odb10/ -m prot -c 24 -sp arabidopsis
cd ..

conda activate genemark

cd predict
mkdir gm3
cd gm3
ln -s ../../round3.maker.output/round3.all.maker.proteins.fasta .
prothint.py --threads 26 --workdir ./g3 ../../data/${genome}.softmask round3.all.maker.proteins.fasta

gmes_petap.pl --EP ./prothint.gff --evidence ./g3/evidence.gff  --cores 26 --max_intron 10000 --seq ../../data/${genome}.softmask  --verbose
cd ..

conda activate maker
mkdir -p snap3
cd snap3

maker2zff -x 0.25 -l 50 -d ../../round3.maker.output/round3_master_datastore_index.log
mv genome.dna rnd3.zff.length50_aed0.25.dna
mv genome.ann rnd3.zff.length50_aed0.25.ann

fathom rnd3.zff.length50_aed0.25.ann rnd3.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom rnd3.zff.length50_aed0.25.ann rnd3.zff.length50_aed0.25.dna -validate > validate.log 2>&1
fathom rnd3.zff.length50_aed0.25.ann rnd3.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..

hmm-assembler.pl rnd3.zff.length50_aed0.25 params > rnd3.zff.length50_aed0.25.hmm

cd ../..

est=
protein=
rmlib=
est2genome=0
protein2genome=0
est_gff=rnd1.all.maker.est2genome.gff
protein_gff=rnd1.all.maker.protein2genome.gff
#rm_gff=rnd1.all.maker.repeats.gff
trna=1
snaphmm=rnd3.zff.length50_aed0.25.hmm
cpus=2
mpiexec=24
always_complete=1
model_org=


maker -OPT
mv maker_opts.ctl rnd4_maker_opts.ctl

sed -i "s/^min_protein=0 /min_protein=30 /g" rnd4_maker_opts.ctl
#sed -i "s/^single_exon=0 /single_exon=1 /g" rnd4_maker_opts.ctl
#sed -i "s/^single_length=250 /single_length=150 /g" rnd4_maker_opts.ctl

sed -i "s/^est_gff= /est_gff=${est_gff} /g" rnd4_maker_opts.ctl
sed -i "s/^protein_gff= /protein_gff=${protein_gff} /g" rnd4_maker_opts.ctl
#####sed -i "s/^trna=0 /trna=${trna} /g" rnd4_maker_opts.ctl
sed -i "s/^snaphmm= /snaphmm=predict\/snap3\/${snaphmm} /g" rnd4_maker_opts.ctl
sed -i "s/^gmhmm= /gmhmm=predict\/gm3\/g3\/GeneMark_ES\/gmhmm.mod /g" rnd4_maker_opts.ctl
sed -i "s/^genome= /genome=data\/${genome} /g" rnd4_maker_opts.ctl
sed -i "s/^est= /est=${est} /g" rnd4_maker_opts.ctl
sed -i "s/^protein= /protein=${protein} /g" rnd4_maker_opts.ctl
sed -i "s/^est2genome=0 /est2genome=${est2genome} /g" rnd4_maker_opts.ctl
sed -i "s/^protein2genome=0 /protein2genome=${protein2genome} /g" rnd4_maker_opts.ctl
sed -i "s/^cpus=1 /cpus=${cpus} /g" rnd4_maker_opts.ctl

sed -i "s/^always_complete=0 /always_complete=${always_complete} /g" rnd4_maker_opts.ctl

sed -i "s/^model_org=all /model_org=${model_org} /g" rnd4_maker_opts.ctl
sed -i "s/^max_dna_len=100000 /max_dna_len=200000 /g" rnd4_maker_opts.ctl
#####sed -i "s/^split_hit=10000 /split_hit=15000 /g" rnd4_maker_opts.ctl

mpiexec -n ${mpiexec} maker -base round4 rnd4_maker_opts.ctl maker_bopts.ctl maker_exe.ctl




#第五轮
cd round4.maker.output
gff3_merge -s -d round4_master_datastore_index.log > rnd4.all.maker.gff
fasta_merge -d round4_master_datastore_index.log
gff3_merge -n -s -d round4_master_datastore_index.log > rnd4.all.maker.noseq.gff
gff3_merge -g -s -n -d round4_master_datastore_index.log > rnd4.gene_model.maker.noseq.gff
cat rnd4.gene_model.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }' >rnd4_genesta.txt
cat round4.all.maker.proteins.fasta |grep ">"|cut -d " " -f 3|tr ":" "\t"|awk '{print $0 }' | awk '{if($2<0.5) sum+=1}END{ print NR, sum / NR }' >>rnd4_genesta.txt

conda activate busco
run_BUSCO.py -i round4.all.maker.proteins.fasta  -o busco_rnd4 -l /home/a2431559261/OLDISK/dataset/embryophyta_odb10/ -m prot -c 24 -sp arabidopsis
cd ..

conda activate genemark

cd predict
mkdir gm4
cd gm4
ln -s ../../round4.maker.output/round4.all.maker.proteins.fasta .
prothint.py --threads 26 --workdir ./g4 ../../data/${genome}.softmask round4.all.maker.proteins.fasta

gmes_petap.pl --EP ./prothint.gff --evidence ./g4/evidence.gff  --cores 26 --max_intron 10000 --seq ../../data/${genome}.softmask  --verbose
cd ..

conda activate maker
mkdir -p snap4
cd snap4

maker2zff -x 0.25 -l 50 -d ../../round4.maker.output/round4_master_datastore_index.log
mv genome.dna rnd4.zff.length50_aed0.25.dna
mv genome.ann rnd4.zff.length50_aed0.25.ann

fathom rnd4.zff.length50_aed0.25.ann rnd4.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom rnd4.zff.length50_aed0.25.ann rnd4.zff.length50_aed0.25.dna -validate > validate.log 2>&1
fathom rnd4.zff.length50_aed0.25.ann rnd4.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..

hmm-assembler.pl rnd4.zff.length50_aed0.25 params > rnd4.zff.length50_aed0.25.hmm

cd ../..

est=
protein=
rmlib=
est2genome=0
protein2genome=0
est_gff=rnd1.all.maker.est2genome.gff
protein_gff=rnd1.all.maker.protein2genome.gff
#rm_gff=rnd1.all.maker.repeats.gff
trna=1
snaphmm=rnd4.zff.length50_aed0.25.hmm
cpus=2
mpiexec=24
always_complete=1
model_org=

maker -OPT
mv maker_opts.ctl rnd5_maker_opts.ctl

sed -i "s/^min_protein=0 /min_protein=30 /g" rnd5_maker_opts.ctl
#sed -i "s/^single_exon=0 /single_exon=1 /g" rnd5_maker_opts.ctl
#sed -i "s/^single_length=250 /single_length=150 /g" rnd5_maker_opts.ctl

sed -i "s/^est_gff= /est_gff=${est_gff} /g" rnd5_maker_opts.ctl
sed -i "s/^protein_gff= /protein_gff=${protein_gff} /g" rnd5_maker_opts.ctl
#####sed -i "s/^trna=0 /trna=${trna} /g" rnd5_maker_opts.ctl
sed -i "s/^snaphmm= /snaphmm=predict\/snap4\/${snaphmm} /g" rnd5_maker_opts.ctl
sed -i "s/^gmhmm= /gmhmm=predict\/gm4\/g4\/GeneMark_ES\/gmhmm.mod /g" rnd5_maker_opts.ctl
sed -i "s/^genome= /genome=data\/${genome} /g" rnd5_maker_opts.ctl
sed -i "s/^est= /est=${est} /g" rnd5_maker_opts.ctl
sed -i "s/^protein= /protein=${protein} /g" rnd5_maker_opts.ctl
sed -i "s/^est2genome=0 /est2genome=${est2genome} /g" rnd5_maker_opts.ctl
sed -i "s/^protein2genome=0 /protein2genome=${protein2genome} /g" rnd5_maker_opts.ctl
sed -i "s/^cpus=1 /cpus=${cpus} /g" rnd5_maker_opts.ctl

sed -i "s/^always_complete=0 /always_complete=${always_complete} /g" rnd5_maker_opts.ctl

sed -i "s/^model_org=all /model_org=${model_org} /g" rnd5_maker_opts.ctl
sed -i "s/^max_dna_len=100000 /max_dna_len=200000 /g" rnd5_maker_opts.ctl
#####sed -i "s/^split_hit=10000 /split_hit=15000 /g" rnd5_maker_opts.ctl

mpiexec -n ${mpiexec} maker -base round5 rnd5_maker_opts.ctl maker_bopts.ctl maker_exe.ctl

cd round5.maker.output
gff3_merge -s -d round5_master_datastore_index.log > rnd5.all.maker.gff
fasta_merge -d round5_master_datastore_index.log
gff3_merge -n -s -d round5_master_datastore_index.log > rnd5.all.maker.noseq.gff
gff3_merge -g -s -n -d round5_master_datastore_index.log > rnd5.gene_model.maker.noseq.gff
cat rnd5.gene_model.maker.noseq.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }' >rnd5_genesta.txt

conda activate busco
run_BUSCO.py -i round5.all.maker.proteins.fasta  -o busco_rnd5 -l /home/a2431559261/OLDISK/dataset/embryophyta_odb10/ -m prot -c 24 -sp arabidopsis
cat round5.all.maker.proteins.fasta |grep ">"|cut -d " " -f 3|tr ":" "\t"|awk '{print $0 }' | awk '{if($2<0.5) sum+=1}END{ print NR, sum / NR }' >>rnd5_genesta.txt

cd ..
