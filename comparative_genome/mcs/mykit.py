#!/usr/bin/python3

import os
import sys
import re
import subprocess
import time
import multiprocessing
import gzip
import pysam
import math
import numpy
import rpy2.robjects as robjects

print('mykit imported')
def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		pass
	try:
		import unicodedata
		unicodedata.numeric(s)
		return True
	except (TypeError, ValueError):
		pass
	return False
	
def multiCMD(cmd,paralist,threads):
	def worker(control,cmd,para):
		control.acquire()
		print('work at {para}'.format(para=para))
		statusoutput=subprocess.getstatusoutput('{cmd} {para}'.format(cmd=cmd,para=para))
		status=statusoutput[0]
		output=statusoutput[1]
		if status != 0:
			with open('workstate','a+') as f:
				f.write('{para}\t{status}\n'.format(para=para,status=status))
			with open('errorlog.log','a+') as f :
				f.write(output)
		if status == 0:
			with open('workstate','a+') as f:
				f.write('{para}\t{status}\n'.format(para=para,status=status))
			with open('log.log','a+') as f:
				f.write(output)
		control.release()
	control=multiprocessing.Semaphore(threads)
	with open('log.log','w') as f:
		pass
	with open('workstate','w') as f:
		pass
	with open('error.log','w') as f:
		pass
	for para in paralist:
		process=multiprocessing.Process(target=worker,args=(control,cmd,para))
		process.start()
	process.join()

def readlist(inputfile):
	lists=[]
	for i in open(inputfile):
		lists.append(i.replace('\n',''))
	return lists

def blastdict(blastresultfile):
	with open(blastresultfile) as f:
		data=f.readlines()
	blastresultdict={}
	for i in data:
		if not blastresultdict.get(i.split('\t')[0],'') :
			blastresultdict[i.split('\t')[0]]=[]
		blastresultdict[i.split('\t')[0]].append(i.split('\t')[1])
	return blastresultdict


def blastdict_incvalue(blastresultfile):
	with open(blastresultfile) as f:
		data=f.readlines()
	blastresultdict={}
	for i in data:
		if not blastresultdict.get(i.split('\t')[0],'') :
			blastresultdict[i.split('\t')[0]]={}
		#blastresultdict[i.split('\t')[0]][i.split('\t')[1]]=[i.split('\t')[2],i.split('\t')[3]]
		blastresultdict[i.split('\t')[0]][i.split('\t')[1]]=int(i.split('\t')[3])
	return blastresultdict
	

def reverse_blastresultdict(raw_blastresultdict):
	reverse_blastdict={}
	for query_gene in raw_blastresultdict:
		for db_result in raw_blastresultdict[query_gene]:
			if db_result not in reverse_blastdict:
				reverse_blastdict[db_result]={}
			reverse_blastdict[db_result][query_gene]=raw_blastresultdict[query_gene][db_result]
	return reverse_blastdict

def wblastdict_incvalue(blastresultdict,outname):
	output = open(outname,'w') 
	for i in blastresultdict:
		for u in blastresultdict[i]:
			writes=i+'\t'+u+'\t'+str(blastresultdict[i][u])+'\n'
			output.write(writes)
	output.close()

def rawannodict(annofile):
	with open(annofile,'r')as f:
		annoinfo=f.readlines()
	annodict={}
#	print(annoinfo)
	for i in annoinfo:
		genename=i.split('\t')[0].replace('\n','')
		anno=i.split('\t')[1].replace('\n','')
		annodict[genename]=anno
	return annodict

def annolist(rawannodict,blastresultdict):
	annolist=[]
#	print('annolist run')
	for i in rawannodict:
#		print(i)
		if i in blastresultdict:
			for u in blastresultdict[i]:
#				print(u+' '+rawannodict[i])
				annolist.append(u+' '+rawannodict[i])
	return annolist

def annodict(rawannodict,blastresultdict):
	annodict={}
#	print('annodict run')
	for i in rawannodict:
#		print(i)
		if i in blastresultdict:
			for u in blastresultdict[i]:
#				print(rawannodict[i])
				annodict[u]=rawannodict[i]
	return annodict

def annodict_inc(rawannodict,reverse_blastresultdict):
	annodict={}
	for i in rawannodict:
		for u in reverse_blastresultdict:
			if i == reverse_blastresultdict[u] :
				annodict[u]=rawannodict[i]
	return annodict
	
def gff_csep(inputfile):
	gffdict={}
	for i in open(inputfile):
		if i.startswith('#'):
			continue
		if '\tmRNA\t' in i or '\ttranscript\t' in i:
			#print(i)
			start=i.split('\t')[3]
			end=i.split('\t')[4]
			chromosome=i.split('\t')[0]
			p_n=i.split('\t')[6]
			for u in i.split('\t')[8].split(';'):
				if 'ID' in u:
					if ':' in u:
						transcript_ID=u.split('=')[1].split(':')[1].replace('\n','')
					else:
						transcript_ID=u.split('=')[1].replace('\n','')
					gffdict[transcript_ID]=[chromosome,start,end,p_n]
	return gffdict

def gtf_csep(inputfile):
	gtfdict={}
	for i in open(inputfile):
		if '\ttranscript\t' in i:
			start=i.split('\t')[3]
			end=i.split('\t')[4]
			chromosome=i.split('\t')[0]
			p_n=i.split('\t')[6]
			for u in i.split('\t')[8].split(';'):
				if 'transcript_id' in u:
					transcript_ID=u.split('"')[1]
					gtfdict[transcript_ID]=[chromosome,start,end,p_n]
	return gtfdict
		
def needgff(genelist,gffdict):
	needgffdict={}
	for i in genelist:
		i=i.split(' ')[0]
		if i in gffdict:
			needgffdict[i]=gffdict[i]
	return needgffdict

def prodict(genelist,gff_csep):
	promoterdict={}
	for i in genelist:
		i=i.split(' ')[0]
		if i in gff_csep:
			p_n=gff_csep[i][3]
			if p_n == '+':
				pstart=str(int(gff_csep[i][1])-2000)
				pend=gff_csep[i][1]
				chromosome=gff_csep[i][0]
				promoterdict[i]=[chromosome,pstart,pend,p_n]
			if p_n == '-':
				pstart=str(int(gff_csep[i][2])+2000)
				pend=gff_csep[i][2]
				chromosome=gff_csep[i][0]
				promoterdict[i]=[chromosome,pend,pstart,p_n]
	return promoterdict

def upstream_dict(genelist,gff_csep,etc_len):
	promoterdict={}
	etc_len=int(etc_len)
	for i in genelist:
		i=i.split(' ')[0]
		if i in gff_csep:
			p_n=gff_csep[i][3]
			if p_n == '+':
				pstart=str(int(gff_csep[i][1])-etc_len)
				pend=gff_csep[i][1]
				chromosome=gff_csep[i][0]
				chromosome=gff_csep[i][0]
				promoterdict[i]=[chromosome,pstart,pend,p_n]
			if p_n == '-':
				pstart=str(int(gff_csep[i][2])+etc_len)
				pend=gff_csep[i][2]
				chromosome=gff_csep[i][0]
				promoterdict[i]=[chromosome,pend,pstart,p_n]
	return promoterdict

def downstream_dict(genelist,gff_csep,etc_len):
	downstream_dict={}
	etc_len=int(etc_len)
	for i in genelist:
		i=i.split(' ')[0]
		if i in gff_csep:
			p_n=gff_csep[i][3]
			if p_n == '+':
				dend=str(int(gff_csep[i][2])+etc_len)
				dstart=gff_csep[i][2]
				chromosome=gff_csep[i][0]
				downstream_dict[i]=[chromosome,dstart,dend,p_n]
			if p_n == '-':
				dend=str(int(gff_csep[i][1])-etc_len)
				dstart=gff_csep[i][1]
				chromosome=gff_csep[i][0]
				downstream_dict[i]=[chromosome,dend,dstart,p_n]
	return downstream_dict


def readfa(inputfile,keepornot):
	###keep or not means whether keep annotation of the seq or not
	fastadict={}
	name=''
	seq=''
	if inputfile.endswith('.gz'):
		for i in gzip.open(inputfile,'r'):
			i=i.decode()
			if '>' in i:
				if name:
					fastadict[name]=seq
					seq=''
				if keepornot == 'yes' :
					name=i.replace('>','').replace('\n','')
				if keepornot == 'no' :
					name=i.split(' ')[0].split('\t')[0].replace('>','').replace('\n','')
			else:
				seq+=i.replace('\n','').strip()
	else:
		for i in open(inputfile):
			if '>' in i:
				if name:
					fastadict[name]=seq
					seq=''
				if keepornot == 'yes' :
					name=i.replace('>','').replace('\n','')
				if keepornot == 'no' :
					name=i.split(' ')[0].split('\t')[0].replace('>','').replace('\n','')
			else:
				seq+=i.replace('\n','').strip()
	fastadict[name]=seq
	return fastadict

def renamefadict(inputfadict,renamedict):
	newdict={}
	for i in inputfadict:
#		print(i)
		if i in renamedict:
			newdict[renamedict[i]]=inputfadict[i]
		else:
			newdict[i]=inputfadict[i]
	return newdict

def reverse(strs):
	return strs[::-1]


def complement(strs):
	if '\n' in strs:
		strs=strs.replace('\n','')

	dict0={'A':'T','G':'C','C':'G','T':'A','a':'t','g':'c','c':'g','t':'a'}
	strs=reverse(strs)
	complementseq=''
	for i in strs:
		complementseq+=dict0[i]
	return complementseq

def etcseq(gff_csep,inputfile,outformat):
#	posstrs=negstrs=''
	outdict={}
	for i in gff_csep:
		posstrs=negstrs=''
		if gff_csep[i][3] == '+':
			start=gff_csep[i][1]
			end=gff_csep[i][2]
			chromosome=gff_csep[i][0]
			posstrs=chromosome+':'+start+'-'+end
			geneseq=subprocess.getoutput('samtools faidx {inputfile} {posstrs}'.format(inputfile=inputfile,posstrs=posstrs))
			#geneseq=pysam.faidx(inputfile,posstrs)
		if gff_csep[i][3] == '-':
			start=gff_csep[i][1]
			end=gff_csep[i][2]
			chromosome=gff_csep[i][0]
			negstrs=chromosome+':'+start+'-'+end
#			print(negstrs)
			geneseq=subprocess.getoutput('samtools faidx --reverse-complement {inputfile} {negstrs}'.format(inputfile=inputfile,negstrs=negstrs))
#			print('1')
#			print(inputfile)
#			print(negstrs)
#			geneseq=pysam.faidx(inputfile,negstrs,'--reverse-complement')
		seq=''
#		print(geneseq)
		for u in geneseq.split('\n'):
			if '>' in u:
				continue
			else:
				seq+=u.replace('\n','').strip()
		outdict[i]=seq
#		print(seq)
#		print(i)
	if outformat == 'dict':
		return outdict


def annofadict(annomsg,annotype,fadict):
	newdict={}
	if annotype == 'list':
		annolist=annomsg
		for i in annolist:
			genename=i.split(' ')[0]
			if genename in fadict:
				newdict[i]=fadict[genename]
	if annotype == 'dict':
		annodict=annomsg
		for i in annodict:
			genename=i
			if i in fadict:
				newdict[annodict[i]+'-{geneID}'.format(geneID=i)]=fadict[genename]
	return newdict


def wfafile(fadict,outputname):
	if outputname.endswith('.gz'):
		output=gzip.open(outputname,'wb')
		for genename in fadict:
			if not fadict[genename]:
				print(genename+' is empty, ignore it...')
				continue
			output.write(b'>'+genename.encode()+b'\n')
			numstrs=0
			for s in fadict[genename]:
				output.write(s.encode())
				numstrs+=1
				if numstrs > 59:
					output.write(b'\n')
					numstrs=0
				output.write(b'\n')
	else:
		output=open(outputname,'w')
		for genename in fadict:
			if not fadict[genename]:
				print(genename+' is empty, ignore it...')
				continue
			output.write('>'+genename+'\n')
			numstrs=0
			for s in fadict[genename]:
				output.write(s)	
				numstrs+=1
				if numstrs > 59:
					output.write('\n')
					numstrs=0
			output.write('\n')
	output.close()


def wgeneinfo(gff_csep,outfile):
	output=open(outfile,'w')
	output.write('#genename\tchromosome\tstart\tend\tp_n\n')
	for i in gff_csep:
		chromosome=gff_csep[i][0]
		start=gff_csep[i][1]
		end=gff_csep[i][2]
		p_n=gff_csep[i][3]
		output.write(i+'\t'+chromosome+'\t'+start+'\t'+end+'\t'+p_n+'\n')

def getstat(cmd,para):
	stat=subprocess.getstatusoutput('{cmd} {para}'.format(cmd=cmd,para=para))[0]
	return stat

def extfile(filename):
	stat=getstat('ls',filename)
	return stat

class evotree:
	def evotree(inputfile,inputtype):
		outdict=inputfile+'.out'
		if inputtype == 'cod-nucl':
			os.system('megacc -a /home/wangjt/Documents/infer_ML_coding-nucleotide.mao -d {inputfile} -o {outdict}'.format(inputfile=inputfile,outdict=outdict))
		if inputtype == 'non-cod-nucl':
			os.system('megacc -a /home/wangjt/Documents/infer_ML_nucleotide.mao -d {inputfile} -o {outdict}'.format(inputfile=inputfile,outdict=outdict))
		if inputtype == 'prot':
			os.system('megacc -a /home/wangjt/Documents/infer_ML_amino_acid.mao -d {inputfile} -o {outdict}'.format(inputfile=inputfile,outdict=outdict))

	def multievo(workdict,threads):
		def buildevotree(s,sysname,inputtype):
			s.acquire()
			inputfile=sysname+'.meg.fas'
			outdir=sysname
			if inputtype == 'cod-nucl':
				os.system('megacc -a /home/wangjt/Documents/infer_ML_coding-nucleotide.mao -d {inputfile} -o {outdir}/'.format(inputfile=inputfile,outdir=outdir))
			if inputtype == 'non-cod-nucl':
				os.system('megacc -a /home/wangjt/Documents/infer_ML_nucleotide.mao -d {inputfile} -o {outdir}/'.format(inputfile=inputfile,outdir=outdir))
			if inputtype == 'prot':
				os.system('megacc -a /home/wangjt/Documents/infer_ML_amino_acid.mao -d {inputfile} -o {outdir}/'.format(inputfile=inputfile,outdir=outdir))
			s.release()
		s=multiprocessing.Semaphore(threads)
		for name in workdict:
			p=multiprocessing.Process(target=buildevotree,args=(s,workdict[name][0],workdict[name][1]))
			p.start()

class cluster:
	def cdhitresdict(inputfile):
		outdict={}
		for i in open(inputfile,'r'):
			if re.search('^>',i):
				name=i.replace('\n','').replace('>','')
				outdict[name]=[]
			else:
				strs=i.split('>')[1].split('.')[0]
				outdict[name].append(strs)
		return outdict
	def sortedcdhitres(inputdict):
		lendict={}
		for i in inputdict:
			lendict[i]=len(inputdict[i])
		sorteddict=dict(sorted(lendict.items(),reverse=True,key=lambda a:a[1]))
		return sorteddict
	def rename_cluster(cdhitresdict,cdhitresfa):
		fadict=readfa(cdhitresfa,'no')
		renamedfadict={}
		for i in cdhitresdict:
			for u in cdhitresdict[i]:
				if u in fadict:
					renamedfadict[i.replace(' ','')]=fadict[u]
		return renamedfadict

class SRAtool:
	def preparedict(inputfile):
		rundict={}
		for i in open(inputfile,'r'):
			i=i.replace('\n','')
			readtype=i.split('\t')[0]
			sysname=i.split('\t')[1]
			if readtype == 'SINGLE':
				rundict[sysname]=[i.split('\t')[0],i.split('\t')[2]]
			if readtype == 'PAIRED':
				rundict[sysname]=[i.split('\t')[0],i.split('\t')[2],i.split('\t')[3]]
		return rundict

	def multiQC(rundict,rawdir):
		def multiQCworker(control,sysname,readtype,rawdir):
			control.acquire()
			if not readtype:
				pass
#				if readtype == 'SINGLE':
#					os.system('java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE {rawdir}/{sysname}.fastq.gz {rawdir}/{sysname}.tri.fastq.gz ILLUMINACLIP:/home/wangjt/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:32'.format(sysname=sysname,rawdir=rawdir))
#					os.system('fastqc -o tmp {rawdir}/{sysname}.fastq.gz {rawdir}/{sysname}.tri.fastq.gz'.format(sysname=sysname,rawdir=rawdir))
#				if readtype == 'PAIRED':
#					os.system('java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE {rawdir}/{sysname}_1.fastq.gz {rawdir}/{sysname}_1.tri.fastq.gz ILLUMINACLIP:/home/wangjt/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:32'.format(sysname=sysname,rawdir=rawdir))
#					os.system('java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE {rawdir}/{sysname}_2.fastq.gz {rawdir}/{sysname}_2.tri.fastq.gz ILLUMINACLIP:/home/wangjt/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:32'.format(sysname=sysname,rawdir=rawdir))
#					os.system('fastqc -o tmp {rawdir}/{sysname}_1.fastq.gz {rawdir}/{sysname}_1.tri.fastq.gz {rawdir}/{sysname}_2.fastq.gz {rawdir}/{sysname}_2.tri.fastq.gz'.format(sysname=sysname,rawdir=rawdir))
			else:
				if readtype == 'SINGLE':
					#os.system('java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE {rawdir}/{sysname}.fastq.gz {rawdir}/{sysname}.tri.fastq.gz MINLEN:32'.format(sysname=sysname,rawdir=rawdir))
					#os.system('fastqc -o tmp {rawdir}/{sysname}.fastq.gz {rawdir}/{sysname}.tri.fastq.gz'.format(sysname=sysname,rawdir=rawdir))
					os.system('fastqc -o tmp {rawdir}/{sysname}.fastq.gz'.format(sysname=sysname,rawdir=rawdir))
				if readtype == 'PAIRED':
					#os.system('java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE {rawdir}/{sysname}_1.fastq.gz {rawdir}/{sysname}_1.tri.fastq.gz MINLEN:32'.format(sysname=sysname,rawdir=rawdir))	
					#os.system('java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE {rawdir}/{sysname}_2.fastq.gz {rawdir}/{sysname}_2.tri.fastq.gz MINLEN:32'.format(sysname=sysname,rawdir=rawdir))
					#os.system('fastqc -o tmp {rawdir}/{sysname}_1.fastq.gz {rawdir}/{sysname}_1.tri.fastq.gz {rawdir}/{sysname}_2.fastq.gz {rawdir}/{sysname}_2.tri.fastq.gz'.format(sysname=sysname,rawdir=rawdir))
					os.system('fastqc -o tmp {rawdir}/{sysname}_1.fastq.gz {rawdir}/{sysname}_2.fastq.gz'.format(sysname=sysname,rawdir=rawdir))
			control.release()
		control=multiprocessing.Semaphore(15)
		for sysname in rundict:
			readtype=rundict[sysname][0]
			process=multiprocessing.Process(target=multiQCworker,args=(control,sysname,readtype,rawdir))
			process.start()
		process.join()



class RNAseq:
	def readtablehead(inputfile):
		f=open(inputfile)
		header=f.readline()
		return header

	def readquantify(inputfile):
		f2=open(inputfile)
		flen=len(f2.readlines())-1
		f2.close()
		f=open(inputfile)
		global header
		header=f.readline()
		quantdict={}
		for i in range(flen):
			msg=f.readline()
			genename=msg.split('\t')[0]
			quantinfo0=msg.split('\t')[1:]
			quantinfo=[]
			for u in quantinfo0 :
				quantinfo.append(u.replace('\n',''))
			quantdict[genename]=quantinfo
		return quantdict
	def need_quantify_table(quantify_dict,annodict):
		need_quantify_dict={}
		for i in quantify_dict:
			if i in annodict:
				need_quantify_dict[i]=quantify_dict[i]
				need_quantify_dict[i].append(annodict[i])
		return need_quantify_dict
	def write_quantify_table(header,quantify_dict,outputname):
		output=open(outputname,'w')
		output.write(header)
		for i in quantify_dict:
			writes=i+'\t'
			for u in quantify_dict[i]:
				writes+=u+'\t'
			writes=re.sub('\t$','\n',writes)
			output.write(writes)
		output.close()

	def featureCounts(outdir):
		bamfilelist=subprocess.getoutput('find {outdir} -iname "alignments.sorted.bam" '.format(outdir=outdir))
		rootdir=subprocess.getoutput('pwd')
		bamfilestr=''
		for i in bamfilelist.split('\n'):
			os.system('ln -s {rootdir}/{i} {filename}'.format(rootdir=rootdir,i=i,filename=i.split('/')[-2]))
			bamfilestr+=i.split('/')[-2]+' '
		print(bamfilestr)
		os.system('featureCounts -o ./featureCounts_quantify.tmp -T 30 -a stringtie_merged.gtf -g gene_id {bamfilestr}'.format(bamfilestr=bamfilestr))
		for i in bamfilelist.split('\n'):
			os.system('rm {filename}'.format(filename=i.split('/')[-2]))
		quantifyresult=open('featureCounts_quantify.result','w')
		for i in open('featureCounts_quantify.tmp'):
			if '#' in i:
				continue
			writes=''
			writes=i.split('\t')[0]+'\t'
			for u in i.split('\t')[6:]:
				writes+=u+'\t'
			writes=re.sub('\t$','',writes)
			quantifyresult.write(writes)
		quantifyresult.close()
		os.system('rm featureCounts_quantify.tmp*')

	def ghostz_nr(inputfile):
		os.system('mkdir tmp')
		def worker(control,inputfile0,i):
			control.acquire()
			print('working at '+str(i))
			os.system('ghostz aln -i {inputfile} -d /home/wangjt/database/annotation/nr/{database} -o tmp/ghostzresult.{i} -q d -t p -a 10'.format(inputfile=inputfile0,database='nr.fa.split.'+str(i)+'.ghostz',i=i))
			control.release()
		control=multiprocessing.Semaphore(3)
		for i in range(0,407):
			process=multiprocessing.Process(target=worker,args=(control,inputfile,i))
			process.start()
		os.system('cat tmp/ghostzresult* > ghostzresult.tmp')
	
	def quantify_check(quantify_file):
		with open (quantify_file) as f:	
			quantify_result0=f.readlines()

		quantify_result={}
		for msg in quantify_result0:
			quantify_result[msg.split('\t')[0]]=msg.split('\t')[1:]

		for unit in quantify_result:
			listlen=len(quantify_result[unit])
			break

		for unit in quantify_result:
			newlistlen=len(quantify_result[unit])
#			print(listlen)
#			print(newlistlen)
			if listlen != newlistlen:
				print('table lens are not the same, check it.')
				sys.exit(1)

		def listmean(inputlist):
			numsum=0.0
			listlen=len(inputlist)
			for i in inputlist:
				numsum+=float(i.replace('\n',''))
			meanvalue=numsum/listlen
			return meanvalue

		def listDx(inputlist):
			meanvalue=listmean(inputlist)
			upone=0.0
			for i in inputlist:
				upone+=(float(i)-meanvalue)**2
			Dx=upone/len(inputlist)
			return Dx

		def listStdev(inputlist):
			return math.sqrt(listDx(inputlist))
		
		valuelist=[10,100,1000]
		for value in valuelist:
			output=open('outputcheck_'+str(value)+'.result','w')
			for unit in quantify_result:
				if unit.startswith('Geneid'):
					writes=unit+'\t'
					for i in quantify_result[unit]:
						writes+=i.replace('\n','')+'\t'
					writes=re.sub('\t$','\n',writes)
					output.write(writes)
					continue
				if listmean(quantify_result[unit]) <= value:
					continue
				writes=unit+'\t'
				for i in quantify_result[unit]:
					writes+=i.replace('\n','')+'\t'
				writes=re.sub('\t$','\n',writes)
				output.write(writes)
			output.close()
	
	def quantify_sort(inputfile):
		def listmean(inputlist):
			numsum=0.0
			listlen=len(inputlist)
			for i in inputlist:
				numsum+=float(i.replace('\n',''))
			meanvalue=numsum/listlen
			return meanvalue

		def listDx(inputlist):
			meanvalue=listmean(inputlist)
			upone=0.0
			for i in inputlist:
				upone+=(float(i)-meanvalue)**2
			Dx=upone/len(inputlist)
			return Dx

		def listStdev(inputlist):
			stdev=math.sqrt(listDx(inputlist))
			return stdev
		
		def Stdev_div_mean (inputlist):
			stdev=listStdev(inputlist)
			div=stdev/listmean(inputlist)
			return div

		quantifydict={}
		out_quantify_dict={}
		for i in open(inputfile):
			quantifydict[i.split('\t')[0]]=i.split('\t')[1:]
		for gene in quantifydict:
			if gene.startswith('Geneid'):
				continue
			out_quantify_dict[gene]=Stdev_div_mean(quantifydict[gene])
		sorted_quantify_dict=dict(sorted(out_quantify_dict.items(),reverse=False,key=lambda a:a[1]))
		return sorted_quantify_dict

	def w_new_sorted_dict(inputfile,outname):
		sortedict=RNAseq.quantify_sort(inputfile)
		output=open(outname,'w')
		with open(inputfile) as f:
			output.write(f.readline())
			resmsg=f.readlines()
			resdict={}
			for i in resmsg:
				resdict[i.split('\t')[0]]=i
		for i in sortedict:
			print(i)
			output.write(resdict[i])
			output.flush()
		output.close()
	
	def single_align(dnafile,transcript_prefix,output_prefix):
		threads=35
		if extfile('tmp'):
			os.system('mkdir tmp')
		if extfile('tmp/{dnafile}.hisat2.index*'.format(threads=threads,dnafile=dnafile)):
			os.system('hisat2-build -p {threads} {dnafile} tmp/{dnafile}.hisat2.index'.format(threads=threads,dnafile=dnafile))
		transcript_1=transcript_prefix+'_R1.fastq.gz'
		transcript_2=transcript_prefix+'_R2.fastq.gz'
		sysname='self_blast_transcript'
		outdir='outfile'
		workdir='workfile'
		
		os.system('run_rnacocktail.py align --align_idx tmp/{dnafile}.hisat2.index --outdir {outdir} --workdir {workdir} --1 {transcript_1} --2 {transcript_2} --threads {threads} --sample {sysname} '.format(dnafile=dnafile,outdir=outdir,workdir=workdir,transcript_1=transcript_1,transcript_2=transcript_2,threads=threads,sysname=sysname))
		os.system('run_rnacocktail.py reconstruct --alignment_bam {outdir}/hisat2/{sysname}/alignments.sorted.bam --outdir {outdir} --workdir {workdir} --threads {threads} --sample {sysname}'.format(outdir=outdir,sysname=sysname,workdir=workdir,threads=threads))

		os.system('gffread {outdir}/stringtie/{sysname}/transcripts.gtf -g {dnafile} -w All.Exon.fa'.format(sysname=sysname,outdir=outdir,dnafile=dnafile))
	
	def stringtie_genemap(input_exon_file):
		genemap={}
		fadict=readfa(input_exon_file,'no')
		for exon_name in fadict:
			#print(exon_name)
			gene_name = exon_name.split('.')[0]+'_'+exon_name.split('.')[1]
			if gene_name not in genemap:
					genemap[gene_name]={}
			genemap[gene_name][exon_name]=len(fadict[exon_name])
#		print(genemap)
		return genemap
		
	def longest_exon(input_exon_file):
		gene_len_map=RNAseq.stringtie_genemap(input_exon_file)
#		print(gene_len_map)
		new_genemap={}
		for gene in gene_len_map:
			longest_exon_len=0
			for exon in gene_len_map[gene]:
					exon_len=gene_len_map[gene][exon]
					if longest_exon_len < exon_len:
							longest_exon_len = exon_len
							new_genemap[gene]=exon
		fadict=readfa(input_exon_file,'no')
		newfadict={}
		for gene in new_genemap:
			newfadict[new_genemap[gene]]=fadict[new_genemap[gene]]
		wfafile(newfadict,input_exon_file+'.longest.fa')
	def cds_genemap(input_cds_file):
		genemap={}
		fadict=readfa(input_cds_file,'no')
		for cds_name in fadict:
			if len(cds_name.split('.')) == 1:
				gene_name=cds_name
			else:
				tailmsg=cds_name.split('.')[-1]
				#print(tailmsg)
			#gene_name=cds_name.split('.')[0]
				#gene_name=cds_name.replace('tailmsg','')
				gene_name=re.sub('\.$','',re.sub(tailmsg+'$','',cds_name))
			if gene_name not in genemap:
				genemap[gene_name]={}
			genemap[gene_name][cds_name]=len(fadict[cds_name])
		#print(genemap)
		return genemap
			
	def longest_cds(input_cds_file):
		gene_len_map=RNAseq.cds_genemap(input_cds_file)
		new_genemap={}
		for gene in gene_len_map:
			longest_cds_len=0
			for cds in gene_len_map[gene]:
					cds_len=gene_len_map[gene][cds]
					if longest_cds_len < cds_len:
							longest_cds_len = cds_len
							new_genemap[gene]=cds
		fadict=readfa(input_cds_file,'no')
		newfadict={}
		for gene in new_genemap:
			newfadict[new_genemap[gene]]=fadict[new_genemap[gene]]
		wfafile(newfadict,input_cds_file+'.longest.fa')
		
class circos:
	def masked_density(input_Repeat_masked,bootstrap):
		density_dict={}
		seqname=''
		file_input=open(input_Repeat_masked)
		while True:
			character = file_input.read(1)
			if character == '':
				break
			if character == '>':
				seqname=file_input.readline().replace('\n','')
				density_dict[seqname]={}
				print(seqname)
				lower_counter=0
				counter=0
				length_counter=0
				continue
			if character.islower():
				lower_counter+=1
			counter+=1
			length_counter+=1
			if counter > int(bootstrap)-1 :
				density_name = str(length_counter-counter)+'-'+str(length_counter)
	#			print(density_name)
				density_dict[seqname][density_name]=lower_counter/counter
				counter=0
				lower_counter=0
		return density_dict

class genome:
	def Repeat_density(input_Repeat_gff):
		seqname_list=[]
		with open(input_Repeat_gff) as f:
			RMgff_msg=f.readlines()
		for i in RMgff_msg:
			if i.startswith('#'):
				continue
			if i.split('\t')[0] not in seqname_list:
				seqname_list.append(i.split('\t')[0])
				with open(i.split('\t')[0]+'.density','w') as f:
					pass
		arange=range(0,0)
		
		for i in RMgff_msg:
			if i.startswith('#'):
				continue
#			if 'Motif:centermere' in i.split('\t')[8]:
			seqname=i.split('\t')[0]
			start=int(i.split('\t')[3])
			end=int(i.split('\t')[4])
			if start not in arange:
				arange=range(start,end)
				with open(seqname+'.density','a+') as f:
					f.write(i.split('\t')[3]+'\n')
		
		if extfile('result'):
			os.system('mkdir result')
#		print(seqname_list)
		for seqname in seqname_list:
#			print(seqname)
			if extfile(seqname+'.density'):
				continue
#			print(seqname)
			with open(seqname+'.density') as f:
				msg=f.readlines()
			if len(msg) <= 3 :
				os.system('rm '+seqname+'.density')
				continue
			print(seqname)
			robjects.r('data<-read.table("{seqname}.density")\nlocation<-data$V1\n png(file="result/{seqname}_centermere_density.png",width=960,height=640)\nplot(density(location))\n dev.off()'.format(seqname=seqname))
#			os.system('rm '+seqname+'.density')
		
	def Repeat_percentage(input_Repeat_masked):
		percentage_dict={}
		seqname=''
		file_input=open(input_Repeat_masked)
		while True:
			character = file_input.read(1)
			if character == '':
				break
			if character == '>':
				if seqname != '':
					percentage_dict[seqname]=low_character/all_character
				seqname=file_input.readline().replace('\n','')
				print(seqname)
				low_character=0
				all_character=0
				continue
			if character.islower():
				low_character+=1
			all_character+=1
		return percentage_dict
	def Repeat_percentage_hardmask(input_Repeat_masked):
		percentage_dict={}
		seqname=''
		file_input=open(input_Repeat_masked)
		while True:
			character = file_input.read(1)
			if character == '':
				break
			if character == '>':
				if seqname != '':
					percentage_dict[seqname]=low_character/all_character
				seqname=file_input.readline().replace('\n','')
				print(seqname)
				low_character=0
				all_character=0
				continue
			if character == 'N':
				low_character+=1
			all_character+=1
		return percentage_dict
			
			
	def genome_transcript_circos(sysname1,sysname2,transcript_prefix,outdir,threads):
		genomefile1=sysname1+'.dna.fa'
		genomefile2=sysname2+'.dna.fa'
		if extfile('tmp'):
			os.system('mkdir tmp')
		syslist=[sysname1,sysname2]
		MCscanX_gff_dict={}
		fadict={}
		arange=range(0,0)
		for sysname in syslist:
			dnafile=sysname+'.dna.fa'
			if extfile('tmp/{dnafile}.1.ht2'.format(dnafile=dnafile)):
				os.system('hisat2-build -p {threads} {dnafile} tmp/{dnafile}'.format(threads=threads,dnafile=dnafile))
			if extfile('tmpworkdir/hisat2/{dnafile}_tmpannotation'.format(dnafile=dnafile)):
				os.system('mkdir -p tmpworkdir/hisat2/{dnafile}_tmpannotation'.format(dnafile=dnafile))
			print('hisat2 --dta --rg-id hisat2 --rg SM:{dnafile}_tmpannotation --threads {threads} -x tmp/{dnafile} -1 {transcript_prefix}_R1.fastq.gz -2 {transcript_prefix}_R2.fastq.gz -S tmpworkdir/hisat2/{dnafile}_tmpannotation/alignments.sam '.format(dnafile=dnafile,threads=threads,transcript_prefix=transcript_prefix))#os.system('run_rnacocktail.py align --align_idx tmp/{dnafile} --outdir {outdir} --workdir tmpworkdir --1 {transcript_prefix}_R1.fastq.gz --2 {transcript_prefix}_R2.fastq.gz --threads {threads} --sample {dnafile}_tmpannotation'.format(dnafile=dnafile,outdir=outdir,transcript_prefix=transcript_prefix,threads=threads))
			#os.system('hisat2 --threads {threads}  -x tmp/{dnafile} -1 {transcript_prefix}_R1.fastq.gz -2 {transcript_prefix}_R2.fastq.gz -S tmpworkdir/hisat2/{dnafile}_tmpannotation/alignments.sam '.format(dnafile=dnafile,threads=threads,transcript_prefix=transcript_prefix))
			os.system('hisat2 --dta --rg-id hisat2 --rg SM:{dnafile}_tmpannotation --threads {threads} -x tmp/{dnafile} -1 {transcript_prefix}_R1.fastq.gz -2 {transcript_prefix}_R2.fastq.gz -S tmpworkdir/hisat2/{dnafile}_tmpannotation/alignments.sam '.format(dnafile=dnafile,threads=threads,transcript_prefix=transcript_prefix))
			os.system('samtools sort -@ {threads} tmpworkdir/hisat2/{dnafile}_tmpannotation/alignments.sam > tmpworkdir/hisat2/{dnafile}_tmpannotation/alignments.sorted.bam'.format(dnafile=dnafile,threads=threads,transcript_prefix=transcript_prefix))
			if extfile('{outdir}/hisat2/{dnafile}_tmpannotation/'.format(outdir=outdir,dnafile=dnafile)):
				os.system('mkdir -p {outdir}/hisat2/{dnafile}_tmpannotation/'.format(outdir=outdir,dnafile=dnafile))
			os.system('cp tmpworkdir/hisat2/{dnafile}_tmpannotation/alignments.sorted.bam {outdir}/hisat2/{dnafile}_tmpannotation/alignments.sorted.bam'.format(outdir=outdir,dnafile=dnafile))
			os.system('stringtie {outdir}/hisat2/{dnafile}_tmpannotation/alignments.sorted.bam  -p {threads} -o {outdir}/stringtie/{dnafile}_tmpannotation//transcripts.gtf -A {outdir}/stringtie/{dnafile}_tmpannotation/gene_abund.tab -v &> tmpworkdir/stringtie.log'.format(outdir=outdir,dnafile=dnafile,threads=threads))
			#os.system('run_rnacocktail.py reconstruct --alignment_bam {outdir}/hisat2/{dnafile}_tmpannotation/alignments.sorted.bam --outdir {outdir} --workdir tmpworkdir --threads {threads} --sample {dnafile}_tmpannotation'.format(outdir=outdir,threads=threads,dnafile=dnafile))
			
			geintro_output=open('{sysname}_geintro.fa'.format(sysname=sysname),'w') 
			used_gene_list=[]
			for gtfmsg in open('{outdir}/stringtie/{dnafile}_tmpannotation/transcripts.gtf'.format(outdir=outdir,threads=threads,dnafile=dnafile)):
				if gtfmsg.startswith('#'):
					continue
				if gtfmsg.split('\t')[2] == 'transcript':
					start=int(gtfmsg.split('\t')[3])
					if start > 2000 :
						start=start-2000
					end=int(gtfmsg.split('\t')[4])+2000
					
					chromosome=gtfmsg.split('\t')[0]
#					if start in arange :
#						continue
#					if end in arange:
#						continue
					real_geneid=gtfmsg.split('\t')[8].split(';')[0].split(' ')[1].replace('"','')
					if real_geneid in used_gene_list:
						continue
					used_gene_list.append(real_geneid)
					arange=range(start,end+1)
					plcmsg=chromosome+':'+str(start)+'-'+str(end)
					GeneID=sysname+'_'+real_geneid+'|'+chromosome+':'+str(start)+'-'+str(end)
					#os.system('samtools faidx {dnafile} {chromosome}:{start}-{end} | sed "s;>.*;>{GeneID};"  >> {dnafile}_geintro.fa'.format(dnafile=dnafile,chromosome=chromosome,start=str(start),end=str(end),GeneID=GeneID))
					seq=pysam.faidx(dnafile,plcmsg).replace(plcmsg,GeneID)
					geintro_output.write(seq)
					MCscanX_gff_dict[GeneID]=[chromosome,str(start),str(end)]

		MCscanX_gff_output=open('MCscanX_{sysname1}_{sysname2}.gff'.format(sysname1=sysname1,sysname2=sysname2),'w')
		for geneID in MCscanX_gff_dict:
			writes=MCscanX_gff_dict[geneID][0]+'\t'+geneID+'\t'+MCscanX_gff_dict[geneID][1]+'\t'+MCscanX_gff_dict[geneID][2]+'\n'
			MCscanX_gff_output.write(writes)
		os.system('makeblastdb -in {sysname1}_geintro.fa -dbtype nucl -parse_seqids -out genomeblast.{sysname1}.db'.format(sysname1=sysname1))
		os.system('blastn -query {sysname2}_geintro.fa -db genomeblast.{sysname1}.db -outfmt 6 -num_threads {threads} -out {sysname1}_{sysname2}.blast.result'.format(sysname1=sysname1,sysname2=sysname2,threads=threads))
		os.system('makeblastdb -in {sysname2}_geintro.fa -dbtype nucl -parse_seqids -out genomeblast.{sysname2}.db'.format(sysname2=sysname2))
		os.system('blastn -query {sysname1}_geintro.fa -db genomeblast.{sysname2}.db -outfmt 6 -num_threads {threads} -out {sysname2}_{sysname1}.blast.result'.format(sysname1=sysname1,sysname2=sysname2,threads=threads))
		
		
	def lastz_axt2MC(axtdir):
	#	filelist=subprocess.getoutput('ls *.axt').split('\n')
		pass
	def axt_inform_extract(inputfile):
		output=open(inputfile+'.circos_link.txt','w')
		for msg in open(inputfile):
			if is_number(msg[0]):
				msglist=msg.split(' ')
				scaffold_A=msglist[1]
				scaffold_A_start=msglist[2]
				scaffold_A_end=msglist[3]
				blast_length=int(scaffold_A_end)-int(scaffold_A_start)
				if blast_length <= 5000:
					continue
#				print(blast_length)
				scaffold_B=msglist[4]
				scaffold_B_start=msglist[5]
				scaffold_B_end=msglist[6]
				writes=scaffold_A+' '+scaffold_A_start+' '+scaffold_A_end+' '+scaffold_B+' '+scaffold_B_start+' '+scaffold_B_end+'\n'
				output.write(writes)
		output.close()
	def multi_axt_inform_extract(axtfile_list,threads):
		def worker_axt_inform_extract(control,inputfile):
			control.acquire()
			print(inputfile)
			genome.axt_inform_extract(inputfile)
			'''
			output=open(inputfile+'.circos_link.txt','w')
			for msg in open(inputfile):
				if is_number(msg[0]):
					msglist=msg.split(' ')
					scaffold_A=msglist[1]
					scaffold_A_start=msglist[2]
					scaffold_A_end=msglist[3]
					scaffold_B=msglist[4]
					scaffold_B_start=msglist[5]
					scaffold_B_end=msglist[6]
					writes=scaffold_A+' '+scaffold_A_start+' '+scaffold_A_end+' '+scaffold_B+' '+scaffold_B_start+' '+scaffold_B_end+'\n'
					output.write(writes)
			output.close()
			'''
			control.release()
			print(inputfile+' has done')
		control=multiprocessing.Semaphore(threads)
		for axtfile in axtfile_list:
#			print(axtfile)
			process=multiprocessing.Process(target=worker_axt_inform_extract,args=(control,axtfile))
			process.start()
		print('check point')
		process.join()
		print('all done')	
			
	def genome_blast_MC(genome_blastresult,outprefix):
		with open(genome_blastresult) as f:
			msglist=f.readlines()
		pesudogfflist=[]
		newmsglist=[]
		sourcechrlist=[]
		aimchrlist=[]
		prefix1='A'
		prefix2='B'
		counter=0
		for msg in msglist:
			msgsplit=msg.split('\t')
			pesudoname1=prefix1+str(counter)
			pesudoname2=prefix2+str(counter)
			pesudo1=msgsplit[0]+'\t'+pesudoname1+'\t'+msgsplit[6]+'\t'+msgsplit[7]+'\n'
			pesudo2=msgsplit[1]+'\t'+pesudoname2+'\t'+msgsplit[8]+'\t'+msgsplit[9]+'\n'
			if msgsplit[0] not in sourcechrlist:
				sourcechrlist.append(msgsplit[0])
			if msgsplit[1] not in aimchrlist:
				aimchrlist.append(msgsplit[1])
			pesudogfflist.append(pesudo1)
			pesudogfflist.append(pesudo2)
			newmsg=msg.replace(msgsplit[0],pesudoname1).replace(msgsplit[1],pesudoname2).replace(msgsplit[6],'1').replace(msgsplit[7],str(int(msgsplit[7])-int(msgsplit[6]))).replace(msgsplit[8],'1').replace(msgsplit[9],str(int(msgsplit[9])-int(msgsplit[8])).replace('-',''))
			newmsglist.append(newmsg)
			counter+=1
		outgff=open(outprefix+'.gff','w')		
		outblast=open(outprefix+'.blast','w')
		for i in pesudogfflist:
			outgff.write(i)
			outgff.flush()
		for i in newmsglist:
			outblast.write(i)
			outblast.flush()
		outgff.close()
		outblast.close()
		os.system('MCScanX '+outprefix+'> /dev/null')
		with open ('{outprefix}.dot.ctl'.format(outprefix=outprefix),'w') as f:
			f.write('2000\n2000\n')
			writes=''
			for i in sourcechrlist:
				writes+=i+','
			writes=re.sub(',$','\n',writes)
			f.write(writes)
			writes=''
			for i in aimchrlist:
				writes+=i+','
			writes=re.sub(',$','\n',writes)
			f.write(writes)
		with open('{outprefix}.circle.ctl'.format(outprefix=outprefix),'w') as f:	
			f.write('2000\n')
			writes=''
			for i in sourcechrlist:
				writes+=i+','
			for i in aimchrlist:
				writes+=i+','
			writes=re.sub(',$','\n',writes)
			f.write(writes)
		pwd=subprocess.getoutput('pwd')
		os.chdir('/home/wangjt/software/MCScanX/downstream_analyses')
		os.system('java dot_plotter -g {pwd}/{outprefix}.gff -s {pwd}/{outprefix}.collinearity -c {pwd}/{outprefix}.dot.ctl -o {pwd}/{outprefix}.dot.PNG'.format(pwd=pwd,outprefix=outprefix))
		os.system('java circle_plotter -g {pwd}/{outprefix}.gff -s {pwd}/{outprefix}.collinearity -c {pwd}/{outprefix}.circle.ctl -o {pwd}/{outprefix}.circle.PNG'.format(pwd=pwd,outprefix=outprefix))
		os.chdir(pwd)
		

	def gtf_trans2genemap(gtffile):
		grpdict={}
		for msg in open(gtffile):
			if '\ttranscript\t' not in msg:
				continue
			for u in msg.split('\t')[8].split(';'):
				if 'gene_id' in u:
					geneid=u.split('"')[1]
				if 'transcript_id' in u:
					transcript_id=u.split('"')[1]
			if geneid not in grpdict:
				grpdict[geneid]=[]
			if transcript_id not in grpdict[geneid]:
				grpdict[geneid].append(transcript_id)
		return grpdict
	
	def bwa_mapping(prefix,indexpath,outdir,preparedict):
		for sysname in preparedict:
			print('work at {prefix} {sysname}'.format(prefix=prefix,sysname=sysname))
			if preparedict[sysname][0] == 'SINGLE':
				ufile=preparedict[sysname][1]
				os.system('bwa mem -t 35 {indexpath} {ufile} > {outdir}/{prefix}.{sysname}.sam'.format(indexpath=indexpath,outdir=outdir,prefix=prefix,sysname=sysname,ufile=ufile))
			if preparedict[sysname][0] == 'PAIRED':
				file1=preparedict[sysname][1]
				file2=preparedict[sysname][2]
				os.system('bwa mem -t 35 {indexpath} {file1} {file2} > {outdir}/{prefix}.{sysname}.sam'.format(indexpath=indexpath,outdir=outdir,prefix=prefix,sysname=sysname,file1=file1,file2=file2))
			os.system('samtools sort -@ 35 {outdir}/{prefix}.{sysname}.sam > {outdir}/{prefix}.{sysname}.bam && samtools index -b -@ 20 {outdir}/{prefix}.{sysname}.bam'.format(outdir=outdir,prefix=prefix,sysname=sysname))
			os.system('rm {outdir}/{prefix}.{sysname}.sam'.format(outdir=outdir,prefix=prefix,sysname=sysname))

	def bowtie2_mapping(prefix,indexpath,outdir,preparedict):
		for sysname in preparedict:
			print('work at {prefix} {sysname}'.format(prefix=prefix,sysname=sysname))
			if preparedict[sysname][0] == 'SINGLE':
				ufile=preparedict[sysname][1]
				os.system('bowtie2 -p 35 -x {indexpath} -U {ufile} -S {outdir}/{prefix}.{sysname}.sam'.format(indexpath=indexpath,outdir=outdir,prefix=prefix,sysname=sysname,ufile=ufile))
			if preparedict[sysname][0] == 'PAIRED':
				file1=preparedict[sysname][1]
				file2=preparedict[sysname][2]
				os.system('bowtie2 -p 35 -x {indexpath} -1 {file1} -2 {file2} -S {outdir}/{prefix}.{sysname}.sam'.format(indexpath=indexpath,outdir=outdir,prefix=prefix,sysname=sysname,file1=file1,file2=file2))
			os.system('samtools sort -@ 35 {outdir}/{prefix}.{sysname}.sam > {outdir}/{prefix}.{sysname}.bam && samtools index -b -@ 20 {outdir}/{prefix}.{sysname}.bam'.format(outdir=outdir,prefix=prefix,sysname=sysname))
			os.system('rm {outdir}/{prefix}.{sysname}.sam'.format(outdir=outdir,prefix=prefix,sysname=sysname))
	
	def mappingrate(filename,seqname):
#		print(filename)
		print(filename)
		bam_reader=pysam.AlignmentFile(filename)
		mappednum=bam_reader.mapped
		unmappednum=bam_reader.unmapped
		sumreads=mappednum+unmappednum
		neednum=bam_reader.count(seqname)
		need_percentage=neednum/sumreads
		return need_percentage
	
	def dictmappingrate(filelist,seqname):
		map_dict={}
		for filen in filelist:
			map_dict[filen]=genome.mappingrate(filen,seqname)
		return map_dict
	
	def mClustermappingrate(bamfile_list,Cluster_list):
		mappingdict={}
		mappingdict['filename']=[]
		for cluster in Cluster_list:
			mappingresult=genome.dictmappingrate(bamfile_list,cluster)
			mappingdict['filename'].append(cluster)
			for filename in mappingresult:
				if filename not in mappingdict:
					mappingdict[filename]=[]
				mappingdict[filename].append(mappingresult[filename])
		return mappingdict
	
	def wMappingdict(mappingdict,outfilename):
		output=open(outfilename,'w')
#		print(mappingdict)
		for i in mappingdict:
			output.write(i+'\t')
			writes=''
			for u in mappingdict[i]:
				writes+=str(u)+'\t'
#			print(writes)
			writes=re.sub('\t$','\n',writes)
			output.write(writes)
			output.flush()

	def multiMpileup (bamlist,refseq,threads):
		def worker(control,bamfile,refseq):
			control.acquire()
			print('work at '+bamfile)
			os.system('samtools mpileup -f {refseq} {bamfile} | gzip > {bamfile}.mpileup.gz'.format(bamfile=bamfile,refseq=refseq))
			control.release()
		control=multiprocessing.Semaphore(threads)
		for bamfile in bamlist:
			process=multiprocessing.Process(target=worker,args=(control,bamfile,refseq))
			process.start()
		process.join()

	def genomeblast(inputfile,aimdb,outname,threads):
		os.system('makeblastdb -in {aimdb} -out tmp -parse_seqids -dbtype nucl'.format(aimdb=aimdb))
		os.system('blastn -query {inputfile} -db tmp -outfmt 6 -out {outname} -num_threads {threads} '.format(inputfile=inputfile,outname=outname,threads=threads))
		os.system('rm tmp*')
	
	def maxgenomeblast(blastresult,length_v,id_v,outfile):
		output=open(outfile,'w')
		os.system("awk '$3>{id_v}' {blastresult}|awk '$4>{length_v}'| sort -n -k 1,1 -k 2,2 -k7,7 > {blastresult}.awk{id_v}_{length_v}.result".format(id_v=id_v,blastresult=blastresult,length_v=length_v))
		utgchrdict={}
		for msg in open('{blastresult}.awk{id_v}_{length_v}.result'.format(id_v=id_v,blastresult=blastresult,length_v=length_v)):
			utgname=msg.split('\t')[0]
			if utgname not in utgchrdict:
				utgchrdict[utgname]={}
			chrname=msg.split('\t')[1]
			if chrname not in utgchrdict[utgname]:
				utgchrdict[utgname][chrname]=[]
			utgchrdict[utgname][chrname].append(msg)
		for utgname in utgchrdict:
			oldblastnum=0
			for chrname in utgchrdict[utgname]:
				newblastnum=len(utgchrdict[utgname][chrname])
				if newblastnum > oldblastnum:
					oldblastnum=newblastnum
					endchr=chrname
			for msg in utgchrdict[utgname][endchr]:
				output.write(msg)
				output.flush()
		os.system('rm {blastresult}.awk{id_v}_{length_v}.result'.format(id_v=id_v,blastresult=blastresult,length_v=length_v))

class genome_annotation:
	def maker_gff_pre(input_gff_file):
		output=open(input_gff_file+'3','w')
		for msg in open(input_gff_file):
			if msg.startswith('>'):
				break
			output.write(msg)
	def gff_complete_reader(input_gff3_file):
		gff_complete_dict={}
		for msg in open(input_gff3_file):
#			print(msg)
			if msg.startswith('#'):
				continue
			msg=msg.split('\t')
			chromosome=msg[0]
			source=msg[1]
			anno_type=msg[2]
			start=msg[3]
			end=msg[4]
			score=msg[5]
			strand=msg[6]
			phase=msg[7]
			raw_attributes=msg[8].split(';')
			attributes={}
			for attributes_msg in raw_attributes:
				attributes[attributes_msg.split('=')[0]]=attributes_msg.split('=')[1].replace('\n','')
			ID=attributes['ID']
			gff_complete_dict[ID]=[chromosome,source,anno_type,start,end,score,strand,phase,attributes]
		return gff_complete_dict
		

class annotation:
	def ghostz_swiss(inputfile):
		os.system('mkdir tmp')
		os.system('ghostz aln -i {inputfile} -d /home/wangjt/database/annotation/swissprot/uniprot_sprot.ghostz -o tmp/swissprot.ghostz.result -q d -t p -a 30'.format(inputfile=inputfile))
		ghostzdict0=blastdict('tmp/swissprot.ghostz.result')
		print(ghostzdict0)
		ghostzdict={}
		for i in ghostzdict0:
			ghostzdict[i]=[]
			for u in range(len(ghostzdict0[i])):
				ghostzdict[i].append(ghostzdict0[i][u].split('|')[1])
#		print(ghostzdict)
		with open('/home/wangjt/database/annotation/swissprot/idmapping.genename.dat','r') as f:
			genenameannodict={}
			while True:
				i=f.readline()
				if not i: break
				genenameannodict[i.split('\t')[0]]=i.split('\t')[2].replace('\n','')
#				print( genenameannodict)
			#print('read dict done')
		needgenenamedict={}
		for i in ghostzdict:
			needgenenamedict[i]=[]
			for u in range(len(ghostzdict[i])):
#				print(ghostzdict[i][u])
				
				if ghostzdict[i][u] in genenameannodict:
					#print('success')
					needgenenamedict[i].append(genenameannodict[ghostzdict[i][u]])
		#print('read done')
		return needgenenamedict

	def ghostz_At(inputfile):
		os.system('mkdir tmp')
		os.system('ghostz aln -i {inputfile} -d /home/wangjt/database/annotation/At/Atanno.ghostz -o tmp/swissprot.ghostz.result -q d -t p -a 30'.format(inputfile=inputfile))
		ghostzdict0=blastdict('tmp/swissprot.ghostz.result')
#		print(ghostzdict0)
		ghostzdict={}
		for i in ghostzdict0:
			ghostzdict[i]=[]
			for u in range(len(ghostzdict0[i])):
				if ghostzdict0[i][u].split('|')[1] in ghostzdict[i]:
					continue
				#print(ghostzdict0[i][u].split('|')[1])
				ghostzdict[i].append(ghostzdict0[i][u].split('|')[1])
		print(ghostzdict)
		return ghostzdict
	def blast_At(inputfile):
		os.system('mkdir tmp')
		os.system('blastn -query {inputfile} -db /home/wangjt/database/annotation/At/Atanno.cds.blast -out tmp/blast.result -outfmt 6 -num_threads 30'.format(inputfile=inputfile))
		blastdict0=blastdict('tmp/blast.result')
		annodict={}
		for msg in open('/home/wangjt/database/rawdb/Brassicease/At/At.cds.annoed.fa'):
			if msg.startswith('>'):
				geneID=msg.split('|')[0].replace('>','').replace(' ','')
				annomsg=msg.split('|')[1]
				annodict[geneID]=annomsg	
		blast_dict={}
		for i in blastdict0:
			blast_dict[i]=[]
			for u in range(len(blastdict0[i])):
				blast_dict[i].append(annodict[blastdict0[i][u]])
		print(blast_dict)
		return blast_dict
	def anno_fadict(ghostzdict,fadict):
		newfadict={}
		for geneID in fadict:
			if geneID in ghostzdict:
				geneID_anno = geneID+' ///'
				for i in ghostzdict[geneID]:
					geneID_anno+=i+' ///'
				geneID_anno=re.sub(' ///$','',geneID_anno)
				newfadict[geneID_anno]=fadict[geneID]
			else:
				newfadict[geneID]=fadict[geneID]
		return newfadict

class transposon:
	def trfallseq(inputfile):
		seqdict={}
		line=0
		for i in open(inputfile):
			line+=1
			if i[0].isdigit() :
				seq1=i.split(' ')[13]
				seq2=i.split(' ')[14].replace('\n','')
				seqdict['trf_'+str(line)+'_seq1']=seq1
				seqdict['trf_'+str(line)+'_seq2']=seq2
		return seqdict
	def trfunitseq(inputfile):
		seqdict={}
		line=0
		for i in open(inputfile):
			line+=1
			if i[0].isdigit():
				seq=i.split(' ')[13]
				seqdict['trfunit_'+str(line)]=seq
		return seqdict
	def readLTRsumresult(inputfile,sysname):
		LTRdict={}
		for i in open(inputfile):
			name=i.split('\t')[0].replace('repeat_region','LTR_retrotransposon')
			annolist=[sysname]
			for u in i.split('\t'):
				u=u.replace('\n','')
				annolist.append(u)
			LTRdict[name]=annolist
		return LTRdict
	def LTR_seq_extractor(input_pass_list,outfile_name,seqfile):
		output=open(outfile_name,'w')
		for msg in open(input_pass_list):
			if msg.startswith('#'):
				continue
			loc_msg=msg.split('\t')[0].replace('..','-')
			seqtype=msg.split('\t')[9]
			seqmsg=pysam.faidx(seqfile,loc_msg)
			output.write(seqmsg)
		output.close()


	'''
	def singleMpileupReader_forStdev(inputfile):
		if inputfile.endswith('.gz'):
			msglist=[]
			for i in gzip.open(inputfile):
				msglist.append(i.decode())
		else:
			with open(inputfile) as f:
				msglist=f.readlines()
		
		outputdict={}
		for i in msglist:
			if i.split('\t')[0] not in outputdict:
				outputdict[i.split('\t')[0]]=[]
			outputdict[i.split('\t')[0]].append(i.split('\t')[3])
		return outputdict
	'''
	def singleMpileupReader_forStdev(bamfile):
		print(extfile(bamfile+'.mpileup.gz'))
		if not extfile(bamfile+'.mpileup'):
			inputfile=bamfile+'.mpileup'
			with open(bamfile+'.mpileup') as f:
				msglist=f.readlines()
		elif not extfile(bamfile+'.mpileup.gz'):
			inputfile=bamfile+'.mpileup.gz'
			msglist=[]
			for i in gzip.open(inputfile):
				msglist.append(i.decode())
		else:
			print('Run samtools mpileup first them run this program!')
#			sys.exit(10)
			time.sleep(10)
		msgdict={}
		lendict={}
		os.system('samtools view -h {bamfile} | grep @SQ | cut -f2,3 > Cluster_len_msg.{bamfile}.tmp'.format(bamfile=bamfile))
		for i in open('Cluster_len_msg.{bamfile}.tmp'.format(bamfile=bamfile)):
			Cluster_name=i.split('\t')[0].split(':')[1]
			lenmsg=i.split('\t')[1].split(':')[1].replace('\n','')
			lendict[Cluster_name]=lenmsg
		os.system('rm Cluster_len_msg.{bamfile}.tmp'.format(bamfile=bamfile))
		outputdict={}
		for i in msglist:
			Cluster_name=i.split('\t')[0]
			if Cluster_name not in outputdict:
				outputdict[Cluster_name]={'len':lendict[Cluster_name],'resmsg':[]}
#				msgdict[Cluster_name]=[]
#				lendict[Cluster_name]=Cluster_len
#			msgdict[Cluster_name].append(i.split('\t')[3])
			outputdict[Cluster_name]['resmsg'].append(i.split('\t')[3])
		return outputdict



	def singleMpileupFilter_forStdev(inputdict):
		def listmean(inputlist):
			numsum=0.0
			listlen=len(inputlist)
			for i in inputlist:
				numsum+=float(i.replace('\n',''))
			meanvalue=numsum/listlen
			return meanvalue
		outdict={}
		for Cluster in inputdict:
			if listmean(inputdict[Cluster]['resmsg']) <= 10:
				continue
			listlen=len(inputdict[Cluster]['resmsg'])
			if listlen/int(inputdict[Cluster]['len']) <= 0.6 :
				continue
			outdict[Cluster]=inputdict[Cluster]['resmsg']
		return outdict
	
	def singleMpileupSort_Stdev(inputdict):
		def listmean(inputlist):
			numsum=0.0
			listlen=len(inputlist)
			for i in inputlist:
				numsum+=float(i.replace('\n',''))
			meanvalue=numsum/listlen
			return meanvalue
		def listDx(inputlist):
			meanvalue=listmean(inputlist)
			upone=0.0
			for i in inputlist:
				upone+=(float(i)-meanvalue)**2
			Dx=upone/len(inputlist)
			return Dx
		def listStdev(inputlist):
			stdev=math.sqrt(listDx(inputlist))
			return stdev
		def Stdev_div_mean (inputlist):
			stdev=listStdev(inputlist)
			div=stdev/listmean(inputlist)
			return div
		out_quantify_dict={}
		for cluster in inputdict:
			out_quantify_dict[cluster]=Stdev_div_mean(inputdict[cluster])
		sorted_quantify_dict=dict(sorted(out_quantify_dict.items(),reverse=False,key=lambda a:a[1]))
		return sorted_quantify_dict	
	
	def mpileup_cover_arrage(input_mpileupfile):
		writes=''
		seq_name=''
		output=open(input_mpileupfile+'.cover.result','w')
		if input_mpileupfile.endswith('.gz'):
			mpileup_read=gzip.open(input_mpileupfile)
		else:
			mpileup_read=open(input_mpileupfile)
		mpileup_read.readline()
		while True:
			if input_mpileupfile.endswith('.gz'):
				msg=mpileup_read.readline().decode()
			else:
				msg=mpileup_read.readline()
			if not msg:
				break
			now_seq_name=msg.split('\t')[0]
			if seq_name != now_seq_name:
				seq_name=now_seq_name
				#!
#				output.write(writes+'\n')
#				writes=''
				output.write('\n>'+seq_name+'\n')
				counter=1
			seq_site=int(msg.split('\t')[1])
			while True:
				if seq_site != counter:
					output.write('0 ')
					#!
#					writes+='0 '
					counter+=1
				else:
					seq_cover=msg.split('\t')[3]
					output.write(seq_cover+' ')
					#!
#					writes+=seq_cover+' '
					counter+=1
					break
#			print(writes)
		output.close()


class blast:
	def uniq_blast_result(input_blastdict):
		unique_gene_list=[]
		for gene in blast_dict:
			if len(blast_dict[gene]) > 1:
				continue
			else:
				unique_gene_list.append(gene)
		return unique_gene_list
	
