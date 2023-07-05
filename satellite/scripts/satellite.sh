#!/bin/bash -l

#-------------------------------------------------------------
# Class-1 Satellite Validation System
# Version 1.0
# Developed by Todd Spindler, 1 April 2019
# Ported to Chen's Hera account with minor fixes, 1 March 2022
#-------------------------------------------------------------

# set up module environment
module purge
module use /scratch2/NCEPDEV/ovp/Lichuan.Chen/modulefiles
module load anaconda-work/1.0.0 mmab/1.0.0
module list

THE_DATE=${1}
STOP_DATE=${2:-$THE_DATE}

TODAY=$(date "+%Y%m%d")
LOG_DATE=${THE_DATE:-$TODAY}

NCORES=12
TASK_QUEUE='batch'
WALL='0:30:00'
PROJ='marine-cpu'
#PROJ='ovp'
LOGPATH=/scratch1/NCEPDEV/stmp2/Samira.Ardani/logs/class-4/satellite
JOB='sat'

mkdir -p $LOGPATH
rm -f $LOGPATH/*.log

export SRC_DIR='/scratch1/NCEPDEV/stmp2/Samira.ardani/github/RTOFS_verif/satellite'

while [[ $THE_DATE -le $STOP_DATE ]]
do
  job1=$(sbatch --parsable -J ${JOB}_1_${THE_DATE} -o $LOGPATH/${JOB}_1_${LOG_DATE}.log -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks=$(($NCORES + 1)) --nodes=1 --wrap "python $SRC_DIR/ush/satellite.py ghrsst $THE_DATE")
  job2=$(sbatch --parsable -J ${JOB}_2_${THE_DATE} -o $LOGPATH/${JOB}_2_${LOG_DATE}.log -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks=$(($NCORES + 1)) --nodes=1 --wrap "python $SRC_DIR/ush/satellite.py aviso $THE_DATE")
  job3=$(sbatch --parsable -J ${JOB}_3_${THE_DATE} -o $LOGPATH/${JOB}_3_${LOG_DATE}.log -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks=$(($NCORES + 1)) --nodes=1 --wrap "python $SRC_DIR/ush/satellite.py smap $THE_DATE")
  job4=$(sbatch --parsable -J ${JOB}_4_${THE_DATE} -o $LOGPATH/${JOB}_4_${LOG_DATE}.log -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks=$(($NCORES + 1)) --nodes=1 --wrap "python $SRC_DIR/ush/satellite.py smos $THE_DATE")
  job5=$(sbatch --parsable --dependency=afterany:${job1}:${job2}:${job3}:${job4} --partition=service -J ${JOB}_trans_${THE_DATE} -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks 1 -o $LOGPATH/transfer_${LOG_DATE}.log --wrap "$SRC_DIR/scripts/satellite_transfer.sh $THE_DATE")
  THE_DATE=$(date --date="$THE_DATE + 1 day" '+%Y%m%d')
done
exit
