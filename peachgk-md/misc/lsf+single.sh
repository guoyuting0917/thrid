#!/bin/sh
#BSUB -q A    # non-project phase1
##BSUB -q YA   # non-project phase2
##BSUB -q PA   # project phase 1
##BSUB -q QA   # project phase 2
#BSUB -W 6000
#BSUB -P [project code]
#BSUB -J single_test

#EXEDIR=`echo $LSF_O_WORKDIR | sed -e 's/^.*\///'`
LOGFILE="log_md"

#cd $LSF_O_WORKDIR

(time ./peachgk_md.out) > $LOGFILE 2>&1

#cd ..
#tar czf $EXEDIR.tar.gz $EXEDIR
