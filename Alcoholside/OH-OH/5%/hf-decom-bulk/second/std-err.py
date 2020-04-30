#!/usr/bin/python3
import pandas as pd 
import numpy as np 
import os
import re
import math


read_dat="thermal.dat"

df=pd.read_csv(read_dat,sep='\s+')

ff=open ('std-err.dat','w')
ff.write('term SD std-err\n')
n=df.shape[0]
sep=5

part=int(n/sep)

i=0
diff2=0
diff=0
for title in list(df):
	i=0
	diff2=0
	diff=0
	while (i<sep):
		a0=df[title].mean()
		df1=df[i*part:(i+1)*part]
		a=df1[title].mean()
		diff=(a-a0)**2
		diff2=diff2+diff
		#print (title,i,a,diff)
		i=i+1

	sd=math.sqrt(diff2/(sep-1))
	se=sd/math.sqrt(sep)
	print (title,sd,se,file=ff)