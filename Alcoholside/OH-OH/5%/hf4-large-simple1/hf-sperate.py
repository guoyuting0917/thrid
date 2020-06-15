#!/usr/bin/python3
import pandas as pd 
import numpy as np 
import os
import re
import math


read_dat="thermal.dat"

df=pd.read_csv(read_dat,sep='\s+')

#ff.write('term SD std-err\n')
n=df.shape[0]
sep=4

part=int(n/sep)

print ('n',n,part)
i=0
diff2=0
diff=0
for title in list(df):
	f=open (title+'.dat','w')
	a1=df[title].mean()
	print (a1,file=f)

	i=0
	diff2=0
	diff=0
	while (i<sep):
		ff=open (title+str(i+1)+'.dat','w') #file name start from 1
		#a0=df.title.mean()
		df1=df[i*part:(i+1)*part]   # from i line to i+1 line, first i have to be 0
		a=df1[title].mean()
		print (a,file=ff)
		#diff2=diff2+diff
		#print (title,i,a)
		i=i+1

	#sd=math.sqrt(diff2/(sep-1))
	#se=sd/math.sqrt(sep)
	#print (title,sd,se,file=ff)