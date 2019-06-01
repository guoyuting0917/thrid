#!/bin/sh

EXEDIR=`pwd | sed -e 's/^.*\///'`
LOGFILE="log_md"
#NPROCS=`wc -l < $MACHINEFILE`

ulimit -Ss 2097152

(time ./peachgk_md.out) > $LOGFILE 2>&1 &

#cd ..
#tar czf $EXEDIR.tar.gz $EXEDIR
