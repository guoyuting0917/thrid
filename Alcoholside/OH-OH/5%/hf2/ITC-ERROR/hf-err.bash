#!/usr/bin/env bash

set -euo pipefail

#a=4

#echo a
#echo ${a}
echo *





for folder in temp_err*
do
	pushd ${folder}
        ./../hf-err.py >hf-err.txt
#	./../plot_fitting.gpi
	#	./../abc1.exe
#	sed -i '1d' data.txt
        popd
done


     
