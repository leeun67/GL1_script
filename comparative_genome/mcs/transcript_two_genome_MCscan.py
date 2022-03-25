#!/usr/bin/python3

import sys
import os
import mykit
import re
import subprocess

sysname1=sys.argv[1]
sysname2=sys.argv[2]
transcript_prefix=sys.argv[3]
outdir=sys.argv[4]
threads=sys.argv[5]

def genome_transcripts_MC(sysname1,sysname2):
	pwd=subprocess.getoutput('pwd')
	os.system("cat {sysname1}_{sysname2}.blast.result | awk '$3>98' | awk '$4>5000' > {sysname1}_{sysname2}.98_5000.blast.result".format(sysname1=sysname1,sysname2=sysname2))
	os.system("cat {sysname2}_{sysname1}.blast.result | awk '$3>98' | awk '$4>5000' > {sysname2}_{sysname1}.98_5000.blast.result".format(sysname1=sysname1,sysname2=sysname2))
	os.system('cp MCscanX_{sysname1}_{sysname2}.gff MCscanX_{sysname2}_{sysname1}.gff'.format(sysname1=sysname1,sysname2=sysname2))
	combine_list=[sysname1+'_'+sysname2,sysname2+'_'+sysname1]
	for combine in combine_list:
		if mykit.extfile(combine):
			os.system('mkdir '+combine)
			os.system('mkdir '+combine+'/query')
			os.system('mkdir '+combine+'/db')
		with open('{combine}.98_5000.blast.result'.format(combine=combine)) as f:
			blastmsg_list=f.readlines()
		query_chromosome_list=[]
		query_msg_dict={}
		db_chromosome_list=[]
		db_msg_dict={}
		for blastmsg in blastmsg_list:
			query_chr=blastmsg.split('|')[1].split(':')[0]
			db_chr=blastmsg.split('|')[2].split(':')[0]
			if query_chr not in query_chromosome_list:
				query_chromosome_list.append(query_chr)
				query_msg_dict[query_chr]=[]
			if db_chr not in db_chromosome_list:
				db_chromosome_list.append(db_chr)
				db_msg_dict[db_chr]=[]
			query_msg_dict[query_chr].append(blastmsg)
			db_msg_dict[db_chr].append(blastmsg)
		for query_chr in query_msg_dict:
			output=open('{combine}/query/{query_chr}.blast'.format(combine=combine,query_chr=query_chr),'w')
			circos_ctl_list=[]
			for msg in query_msg_dict[query_chr]:
				output.write(msg)
			output.close()
			os.system('cp MCscanX_{sysname1}_{sysname2}.gff {combine}/query/{query_chr}.gff'.format(sysname1=sysname1,sysname2=sysname2,combine=combine,query_chr=query_chr))
			os.system('MCScanX {combine}/query/{query_chr} &> /dev/null'.format(combine=combine,query_chr=query_chr))
			for msg in open('{combine}/query/{query_chr}.collinearity'.format(combine=combine,query_chr=query_chr)):
				if msg.startswith('#'):
					continue
				chr1=msg.split('|')[1].split(':')[0]
				chr2=msg.split('|')[2].split(':')[0]
				if chr1 not in circos_ctl_list:
					circos_ctl_list.append(chr1)
				if chr2 not in circos_ctl_list:
					circos_ctl_list.append(chr2)
			with open('{combine}/query/{query_chr}.dot.ctl'.format(combine=combine,query_chr=query_chr),'w') as f:
				f.write('2000\n')
				writes=''
				for nchr in circos_ctl_list:
					writes+=nchr+','
				writes=re.sub(',$','',writes)
				f.write(writes)
			os.chdir('/home/a2431559261/soft/MCScanX/downstream_analyses')
			os.system('java circle_plotter -g {pwd}/MCscanX_{combine}.gff -s {pwd}/{combine}/query/{query_chr}.collinearity -c {pwd}/{combine}/query/{query_chr}.dot.ctl -o {pwd}/{combine}/query/{query_chr}.dot.PNG &> /dev/null'.format(pwd=pwd,combine=combine,query_chr=query_chr))
			os.chdir(pwd)
		for db_chr in db_msg_dict:
			output=open('{combine}/db/{db_chr}.blast'.format(combine=combine,db_chr=db_chr),'w')
			circos_ctl_list=[]
			for msg in db_msg_dict[db_chr]:
				output.write(msg)
			output.close()
			os.system('cp MCscanX_{sysname1}_{sysname2}.gff {combine}/db/{db_chr}.gff'.format(sysname1=sysname1,sysname2=sysname2,db_chr=db_chr,combine=combine))
			os.system('MCScanX {combine}/db/{db_chr} &> /dev/null'.format(combine=combine,db_chr=db_chr))
			for msg in open('{combine}/db/{db_chr}.collinearity'.format(combine=combine,db_chr=db_chr)):
				if msg.startswith('#'):
					continue
				chr1=msg.split('|')[1].split(':')[0]
				chr2=msg.split('|')[2].split(':')[0]
				if chr1 not in circos_ctl_list:
					circos_ctl_list.append(chr1)
				if chr2 not in circos_ctl_list:
					circos_ctl_list.append(chr2)
			with open('{combine}/db/{db_chr}.dot.ctl'.format(combine=combine,db_chr=db_chr),'w') as f:
				f.write('2000\n')
				writes=''
				for nchr in circos_ctl_list:
					writes+=nchr+','
				writes=re.sub(',$','',writes)
				f.write(writes)
			os.chdir('/home/a2431559261/OLDISK/soft/MCScanX/downstream_analyses')
			os.system('java circle_plotter -g {pwd}/MCscanX_{combine}.gff -s {pwd}/{combine}/db/{db_chr}.collinearity -c {pwd}/{combine}/db/{db_chr}.dot.ctl -o {pwd}/{combine}/db/{db_chr}.dot.PNG &> /dev/null'.format(pwd=pwd,combine=combine,db_chr=db_chr))
			os.chdir(pwd)

mykit.genome.genome_transcript_circos(sysname1,sysname2,transcript_prefix,outdir,threads)
genome_transcripts_MC(sysname1,sysname2)
os.system('rm tmp genomeblast*.db* *.blast.result *._geintro.fa tmpworkdir MCscanX_* outfile -rf')
