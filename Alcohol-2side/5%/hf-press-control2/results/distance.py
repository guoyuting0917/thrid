#!/usr/bin/python3
import pandas as pd 
import numpy as np 
import os
import re



input_dat="log.distance"
out_dat="log.distance1"

def thermal(input_dat):
	f=open(out_dat,'w')

	for line in open(input_dat,'r'):

			if "#" in line:
				continue
			else:      
				f.write(line)

	return()	


read_dat=thermal(input_dat)


#pd.set_option('display.max_columns', None)

df=pd.read_csv(out_dat,sep='\s+')
#df['nve']=df['TotEng']+df['f_nvt']
print(df.tail(1000000).mean())
#print(df.describe(include='all'))