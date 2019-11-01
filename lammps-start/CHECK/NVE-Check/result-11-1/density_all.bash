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
        ./../cal.out
#	./../abc.exe
	#	./../abc1.exe
#	sed -i '1d' data.txt
        popd
done

sed -i '1d' $data_file
     
