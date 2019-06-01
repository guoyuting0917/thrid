#!/bin/sh

EXEDIR=`pwd | sed -e 's/^.*\///'`
MPIRUN="/opt/openmpi/bin/mpirun"
LOGFILE="log_md"
OPMPIOPT="-bysocket -bind-to-core"
MACHINEFILE="/opt/tools/openmpihosts"
NPROCS=8

ulimit -Ss unlimited

#(time $MPIRUN -np $NPROCS -hostfile $MACHINEFILE $OPMPIOPT ./peachgk_md.out) > $LOGFILE 2>&1
(time $MPIRUN -np $NPROCS -hostfile $MACHINEFILE $OPMPIOPT ./peachgk_md.out) > $LOGFILE 2>&1 &

#cd ..
#tar czf $EXEDIR.tar.gz $EXEDIR
