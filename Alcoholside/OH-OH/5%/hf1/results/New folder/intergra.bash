#!/usr/bin/env bash

set -euo pipefail

#a=4

#echo a
#echo ${a}
echo *




for folder in bigsys*
do
	pushd ${folder}
        ./../intergra.py
	cat intergra_* > all_intergra.dat
	rm -f intergra_*
#	./../abc.exe
	#	./../abc1.exe
#	sed -i '1d' data.txt
        popd
done


     
