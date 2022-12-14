#!/bin/bash -l

#-------------------------------------------------------
# Florida Cable Transport Class-4 Validation System
# Version 1.1
# Todd Spindler
# 5 June 2018
#-------------------------------------------------------

# set up module environment
module purge
module use /scratch2/NCEPDEV/ovp/Lichuan.Chen/modulefiles
module load anaconda-work/1.0.0 mmab/1.0.0
module list
which python

THE_DATE=${1}
NCORES=12
TASK_QUEUE='batch'
WALL='3:00:00'
PROJ='marine-cpu'
#PROJ='ovp'
LOGPATH=/scratch2/NCEPDEV/stmp1/Lichuan.Chen/logs/class-4/cable
JOB="cable"

mkdir -p $LOGPATH
rm -f $LOGPATH/*.log

export SRCDIR='/scratch2/NCEPDEV/ovp/Lichuan.Chen/VPPPG/Global_RTOFS/EMC_ocean-verification/cable'

#bsub -J ${JOB} -q $TASK_QUEUE -P $PROJ -W $WALL -n $NCORES -R $RUSAGE -x -oo $LOGPATH/${JOB}_${THE_DATE}.log python $CABLE_SRC_DIR/ush/cable_transport.py $THE_DATE

job1=$(sbatch --parsable -J ${JOB}_${THE_DATE} -o $LOGPATH/${JOB}_${THE_DATE}.log -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks=$(($NCORES + 1)) --nodes=1 --wrap "python $SRCDIR/ush/cable_transport.py $THE_DATE")
job2=$(sbatch --parsable --dependency=afterok:$job1 --partition=service -J ${JOB}_transfer_${THE_DATE} -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks 1 -o $LOGPATH/transfer_${THE_DATE}.log --wrap "$SRCDIR/scripts/cable_transfer.sh")

exit
