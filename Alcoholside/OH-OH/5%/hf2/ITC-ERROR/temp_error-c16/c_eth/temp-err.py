#!/usr/bin/python3
import pandas as pd 
import numpy as np 
import scipy.constants as C
import sys
import math
#import matplotlib.pyplot as plt
#from scipy import optimize




def temperr(keyword,part):


	result_dat1="2_temp_weight_"+keyword+"_ave.dat"
	df1=pd.read_csv(result_dat1,sep='\s+')
	df1['temp_tot']=df1['temp']
	i=2
	while (i<=n):
		result_dat2=str(i)+"_temp_weight_"+keyword+"_ave.dat"
		df2=pd.read_csv(result_dat2,sep='\s+')
		df1['temp_tot']=df1['temp_tot']+df2['temp']
		i=i+1

	df1['temp_tot']=df1['temp_tot']/part #n
	print(df1['temp_tot'])

	
	i=2
	df1['sigma']=0
	while (i<=n):
		result_dat2=str(i)+"_temp_weight_"+keyword+"_ave.dat"
		df3=pd.read_csv(result_dat2,sep='\s+')
		df3['temp_err']=df3['temp']-df1['temp_tot']
		df3['temp_err']=df3['temp_err']**2
		df1['sigma']=df1['sigma']+df3['temp_err']
		print(i,df3['temp_err'])
		i=i+1
	#df1['sigma']=df1['sigma']/part
	df1['sigma']=np.sqrt(df1['sigma']/part)


	err_dat=open("err_dat",'w')
	print(df1['sigma'],file=err_dat)


	return()


n=4
#temperr('sio2',n)
temperr('c_eth',n)
