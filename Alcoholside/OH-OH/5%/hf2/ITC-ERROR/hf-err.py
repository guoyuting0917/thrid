#!/usr/bin/python3
import pandas as pd 
import numpy as np 
import scipy.constants as C
import sys
import math
#import matplotlib.pyplot as plt
#from scipy import optimize





def tempjump(part):

	fo=open('hf'+part+'dat','r')
	line=fo.readline().split()
	hf=float(line[3]) #MW/m2
	print('heatflux=',hf)
	itc1=hf
	itc2=hf


	return(itc1,itc2)


n=4
a1,a11=tempjump('1_')
a2,a22=tempjump('2_')
a3,a33=tempjump('3_')
a4,a44=tempjump('4_')

a=(a1+a2+a3+a4)/n
aa=(a11+a22+a33+a44)/n

sigma=math.sqrt(((a1-a)**2+(a2-a)**2+(a3-a)**2+(a4-a)**2)/(n-1))
s=sigma/math.sqrt(n)

sigma1=math.sqrt(((a11-aa)**2+(a22-aa)**2+(a33-aa)**2+(a44-aa)**2)/(n-1))
s1=sigma1/math.sqrt(n)

print('starderror_left=',sigma,'starderror of mean=',s)

print('starderror_right=',sigma1,'starderror of mean=',s1)

