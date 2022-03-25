import numpy as np
import pandas as pd
from scipy import stats,integrate
import matplotlib.pyplot as plt
import seaborn as sns
import sys


insertime='inser_time'
title=sys.argv[1]


tmy=[]
for i in range(0,6000000,100000):
        i=i/1000000.0
        tmy.append(i)





f=open(insertime,'r')
f=f.read()
f=f.split('\n')
del f[-1]



ti_num=[]
for i in f:
        ti_num.append(int(i)/1000000.0)


plt.ylabel("Number of LTRs")
plt.title(title)
plt.grid(False)



sns.distplot(ti_num,bins=tmy,kde=False, rug=False,axlabel='Mya')

plt.savefig('inser.png')
