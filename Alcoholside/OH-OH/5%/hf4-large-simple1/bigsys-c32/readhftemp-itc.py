#!/usr/bin/python3
import pandas as pd 
import numpy as np 
import os
import re
import math

for line in open('ITC.dat'):
	if 'temp' in line:
		btempave=float(line.split()[1])
		ttempave=float(line.split()[4])

f=open('v_hf_bsolv.dat')
hf_bsolvave=f.read()
f=open('v_hf_bsurf.dat')
hf_bsurfave=f.read()

f=open('v_hf_tsolv.dat')
hf_tsolvave=f.read()
f=open('v_hf_tsurf.dat')
hf_tsurfave=f.read()

itc_bsolvave=float(hf_bsolvave)/btempave/1e6
itc_tsolvave=float(hf_tsolvave)/ttempave/1e6
itc_bsurfave=float(hf_bsurfave)/btempave/1e6
itc_tsurfave=float(hf_tsurfave)/ttempave/1e6

print('ave',itc_bsurfave)

diff_bsolv2=0
diff_tsolv2=0
diff_bsurf2=0
diff_tsurf2=0 

sep=4
i=1

while (i<5):
		
	for line in open(str(i)+'_ITC.dat'):

		f=open('v_hf_bsolv'+str(i)+'.dat')
		hf_bsolv=f.read()
		f=open('v_hf_tsolv'+str(i)+'.dat')
		hf_tsolv=f.read()
		f=open('v_hf_bsurf'+str(i)+'.dat')
		hf_bsurf=f.read()
		f=open('v_hf_tsurf'+str(i)+'.dat')
		hf_tsurf=f.read()

		if 'temp' in line:

			btemp=float(line.split()[1])
			ttemp=float(line.split()[4])
			itc_bsolv=float(hf_bsolv)/btemp/1e6
			itc_tsolv=float(hf_tsolv)/ttemp/1e6
			itc_bsurf=float(hf_bsurf)/btemp/1e6
			itc_tsurf=float(hf_bsurf)/ttemp/1e6

			print (i,itc_bsurf)

			diff_bsolv=(itc_bsolv - itc_bsolvave )**2
			diff_bsolv2=diff_bsolv2+diff_bsolv
			
			diff_tsolv=(itc_tsolv - itc_tsolvave )**2
			diff_tsolv2=diff_tsolv2+diff_tsolv

			diff_bsurf=(itc_bsurf - itc_bsurfave )**2
			diff_bsurf2=diff_bsurf2+diff_bsurf
			diff_tsurf=(itc_tsurf - itc_tsurfave )**2
			diff_tsurf2=diff_tsurf2+diff_tsurf

	i=i+1

sd_bsolv=math.sqrt(diff_bsolv2/(sep-1))/math.sqrt(sep)
sd_tsolv=math.sqrt(diff_tsolv2/(sep-1))/math.sqrt(sep)
sd_bsurf=math.sqrt(diff_bsurf2/(sep-1))/math.sqrt(sep)
sd_tsurf=math.sqrt(diff_tsurf2/(sep-1))/math.sqrt(sep)

fff=open('itc_err.dat','w')
fff.write('sd_bsolv sd_tsolv sd_tsurf sd_bsurf\n')

print (sd_bsolv,sd_tsolv,sd_tsurf,sd_bsurf,file=fff)