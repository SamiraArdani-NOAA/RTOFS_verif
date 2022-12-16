#!/bin/bash -l

#-------------------------------------------------------------
# WBC Frontal Analysis Validation System
# Version 1.0
# Developed by Todd Spindler, 24 Aug 2020
# Ported to Chen's Hera account with minor fixes, 1 March 2022
#-------------------------------------------------------------

# set up module environment
module purge
module use /scratch2/NCEPDEV/ovp/Lichuan.Chen/modulefiles
module load anaconda-work/1.0.0 mmab/1.0.0
module list

THE_DATE=${1:-`date --date yesterday +%Y%m%d`}
NCORES=12
TASK_QUEUE='batch'
WALL='0:30:00'
PROJ='ovp'
#PROJ='marine-cpu'
LOGPATH=/scratch2/NCEPDEV/stmp1/Lichuan.Chen/logs/class-4/fronts
JOB="fronts"

mkdir -p $LOGPATH
rm -f $LOGPATH/*.log

export SRCDIR='/scratch2/NCEPDEV/ovp/Lichuan.Chen/VPPPG/Global_RTOFS/EMC_ocean-verification/fronts'

job1=$(sbatch --parsable -J ${JOB}_${THE_DATE} -o $LOGPATH/${JOB}_${THE_DATE}.log -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks=$(($NCORES + 1)) --nodes=1 --wrap "python $SRCDIR/ush/global_fronts.py $THE_DATE")
job2=$(sbatch --parsable --dependency=afterok:$job1 --partition=service -J ${JOB}_transfer_${THE_DATE} -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks 1 -o $LOGPATH/transfer_${THE_DATE}.log --wrap "$SRCDIR/scripts/global_fronts_transfer.sh $THE_DATE")

exit
