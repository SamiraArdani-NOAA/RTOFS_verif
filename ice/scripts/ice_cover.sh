#!/bin/bash -l

#-------------------------------------------------------------
# Global RTOFS Ice Concentration Metrics
# Version 1.0
# Developed by Todd Spindler, 24 Aug 2020
# Ported to Chen's Hera account with minor fixes, 1 March 2022
#-------------------------------------------------------------

# set up module environment
module purge
module use /scratch1/NCEPDEV/stmp2/Samira.Ardani/modulefiles
module load anaconda-work/1.0.0 mmab/1.0.0
module list

if [[ -n $1 ]]; then
  THE_DATE=$1
else
  THE_DATE=`date --date='yesterday' +%Y%m%d`
fi

NCORES=9
TASK_QUEUE='batch'
WALL='0:30:00'
PROJ='ovp'
#PROJ='marine-cpu'
LOGPATH=/scratch1/NCEPDEV/stmp2/Samira.Ardani/logs/ice
JOB="ice"

mkdir -p $LOGPATH
rm -f $LOGPATH/*.log

export SRCDIR='/scratch1/NCEPDEV/stmp2/Samira.Ardani/github/RTOFS_verif/ice'

job1=$(sbatch --parsable -J ${JOB}_${THE_DATE} -o $LOGPATH/${JOB}_${THE_DATE}.log -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks=$(($NCORES + 1)) --nodes=1 --wrap "python $SRCDIR/ush/ice_cover.py $THE_DATE")
job2=$(sbatch --parsable --dependency=afterok:$job1 --partition=service -J ${JOB}_transfer_${THE_DATE} -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks 1 -o $LOGPATH/transfer_${THE_DATE}.log --wrap "$SRCDIR/scripts/ice_cover_transfer.sh $THE_DATE")

exit
