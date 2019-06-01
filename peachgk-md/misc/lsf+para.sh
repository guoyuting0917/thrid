#!/bin/sh
#BSUB -q A    # non-project phase1
##BSUB -q YA   # non-project phase2
##BSUB -q PA   # project phase 1
##BSUB -q QA   # project phase 2
#BSUB -W 6000
##BSUB -P [project code]
#BSUB -J mpi_test
#BSUB -n 16

#EXEDIR=`echo $LSF_O_WORKDIR | sed -e 's/^.*\///'`
LOGFILE="log_md"
NPROCS=16

#cd $LSF_O_WORKDIR

ulimit -Ss unlimited
(time mpirun -v -np $NPROCS dplace -s1 ./peachgk_md.out) > $LOGFILE 2>&1

#cd ..
#tar czf $EXEDIR.tar.gz $EXEDIR
