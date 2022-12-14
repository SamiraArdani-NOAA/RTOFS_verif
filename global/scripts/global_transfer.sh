#!/bin/bash -l
#
# this is the global image transfer task
#

FORMAT=${1:-nc}
FCST=${2:-n048}
THE_DATE=${3:-`date +%Y%m%d`}
PARA=${4:-False}

if [[ $PARA = 'PARA' ]]; then
  LOCAL_DIR="/scratch2/NCEPDEV/stmp1/Lichuan.Chen/para/images/${THE_DATE}"
  REMOTE_DIR="/home/www/polar/develop/global/para/${FORMAT}"
else
  LOCAL_DIR="/scratch2/NCEPDEV/stmp1/Lichuan.Chen/images/mpi/${THE_DATE}"
  REMOTE_DIR="/home/www/polar/global/${FORMAT}"
fi

REMOTE_ID='emc.rtofs@emcrzdm.ncep.noaa.gov'

echo "Starting $FORMAT image transfer at `date`"

for size in large small; do
  ssh ${REMOTE_ID} mkdir -p $REMOTE_DIR/images/$size
  ssh ${REMOTE_ID} mkdir -p $REMOTE_DIR/archive/images/$THE_DATE/$size
  ssh ${REMOTE_ID} rm -f $REMOTE_DIR/images/$size/*${FCST}*.png
  scp -q -o "BatchMode yes" $LOCAL_DIR/$size/*${FCST}*.png ${REMOTE_ID}:$REMOTE_DIR/images/$size/. &
  scp -q -o "BatchMode yes" $LOCAL_DIR/$size/*${FCST}*.png ${REMOTE_ID}:$REMOTE_DIR/archive/images/$THE_DATE/$size/. &
done

wait

echo "Completed $FORMAT image transfer at `date`"

# cleanup
rm -f $LOCAL_DIR/large/*${FCST}*.png
rm -f $LOCAL_DIR/small/*${FCST}*.png

exit
