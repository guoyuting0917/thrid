#!/usr/bin/env bash

set -euo pipefail

#a=4

#echo a
#echo ${a}
echo *

#exit
data_file=data.txt


for folder in bigsys*
do
    
    pushd ${folder}
    if [ ! -d "temp_error" ]; then
	mkdir temp_error
	cp temp*.dat temp_error
    fi
    cp  thermo.log temp_error
    cp  wall.txt temp_error
    cp  region_liquid.dat temp_error
    cp  region_sio2.dat   temp_error


#	./../abc.exe
	#	./../abc1.exe
#	sed -i '1d' data.txt
        popd
done

     
