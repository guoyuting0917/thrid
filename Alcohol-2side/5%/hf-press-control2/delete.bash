#!/usr/bin/env bash

set -euo pipefail

#a=4

#echo a
#echo ${a}


for folder in c*
do pushd ${folder}
   pwd
 #  rm *.gpi
 #  rm *.png
 #  rm *.dat
 #  rm npt05fs*
 #  rm cal_temp*
   rm -f *.dcd
   rm -f log*
   rm -f *log
   rm -f density*
   rm -f in.sys-de*
   rm -f *.dat
   find . -name "*" -type f -size 0c | xargs -n 1 rm -f
   popd
done



	
