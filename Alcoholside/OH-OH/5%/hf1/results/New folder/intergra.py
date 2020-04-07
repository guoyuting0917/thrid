#!/usr/bin/python3
import pandas as pd
import scipy.constants as C
import sys
import numpy as np

def intergration(keyword,a,b,comment):
	read_dat="density_"+keyword+"_ave.dat"
	result_dat=open("intergra_"+keyword+comment+".dat",'w')


	df=pd.read_csv(read_dat,sep='\s+')
	intergra_density_1=df.density[(df['Coord1']>a)&(df['Coord1']<b)].sum()
	print(keyword+'_'+comment,'=',intergra_density_1,file=result_dat)
	return()

	#print(df.temp[(df['Coord1']<b)&(df['Coord1']>a)])	
	#intergra_density_2=df.density[(df['Coord1']>c)&(df['Coord1']<d)].sum()
	#intergra_density_ave=(intergra_density_1+intergra_density_2)/2
	#Ncount_ave=df.Ncount[(df['Coord1']>a)&(df['Coord1']<b)].sum()
	#temp_ave=Ncotemp_ave/Ncount_ave
    #print(df.temp[(df['Coord1']>a)&(df['Coord1']<b)].mean()) #simple average temp 
	
def ave_density(keyword,a,b,comment):
	read_dat="density_"+keyword+"_ave.dat"
	result_dat=open("intergra_"+keyword+comment+".dat",'w')
	df=pd.read_csv(read_dat,sep='\s+')
	den_ave=df.density[(df['Coord1']>a)&(df['Coord1']<b)].mean()
	print(keyword+'_'+comment,'=',den_ave,file=result_dat)
	return()


read_dat="wall.txt"
df=pd.read_csv(read_dat,sep='\s+')
df.wall_l[0]=0

print(df.wall_l[0],df.wall_r[0])

intergration("eth",a=df.wall_l[0]+26,b=df.wall_l[0]+32,comment='wall_l')
intergration("eth",a=df.wall_r[0]-32,b=df.wall_r[0]-26,comment='wall_r')
ave_density("eth",a=df.wall_l[0]+55,b=df.wall_r[0]-55,comment='ave')

intergration("C24",a=df.wall_l[0]+26,b=df.wall_l[0]+32,comment='wall_l')
intergration("C24",a=df.wall_r[0]-32,b=df.wall_r[0]-26,comment='wall_r')
ave_density("C24",a=df.wall_l[0]+55,b=df.wall_r[0]-55,comment='ave')

intergration("h_eth",a=df.wall_l[0]+26,b=df.wall_l[0]+32,comment='wall_l')
intergration("h_eth",a=df.wall_r[0]-32,b=df.wall_r[0]-26,comment='wall_r')

intergration("o_eth",a=df.wall_l[0]+26,b=df.wall_l[0]+32,comment='wall_l')
intergration("o_eth",a=df.wall_r[0]-32,b=df.wall_r[0]-26,comment='wall_r')

intergration("c_eth",a=df.wall_l[0]+26,b=df.wall_l[0]+32,comment='wall_l')
intergration("c_eth",a=df.wall_r[0]-32,b=df.wall_r[0]-26,comment='wall_r')


intergration("C24",a=df.wall_l[0]+32,b=df.wall_l[0]+36.3,comment='wall_l2')
intergration("C24",a=df.wall_r[0]-36.3,b=df.wall_r[0]-32,comment='wall_r2')

intergration("c_eth",a=df.wall_l[0]+32,b=df.wall_l[0]+36.3,comment='wall_l2')
intergration("c_eth",a=df.wall_r[0]-36.3,b=df.wall_r[0]-32,comment='wall_r2')




