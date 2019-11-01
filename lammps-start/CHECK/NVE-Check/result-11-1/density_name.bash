#!/usr/bin/env bash

set -euo pipefail

#a=4

#echo a
#echo ${a}
echo *

#exit
train_file=filename.txt
ave_file=fileave.txt
data_file=data.txt

path=./
for folder in bigsys*
do
    pushd ${folder}
    find $path -name '*.dat' > $train_file
    sed -i 's:./::g' $train_file
    
    find $path -name '*.dat' > $ave_file
    sed -i 's:./::g' $ave_file    
    sed -i 's/.dat/_ave.dat/g' $ave_file
    paste $train_file $ave_file > $data_file
    rm $train_file $ave_file

    
#        ./../cal_den.out
#	./../abc.exe
#	./../abc1.exe
        popd
	done
     
