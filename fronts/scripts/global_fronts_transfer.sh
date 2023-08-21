#!/bin/bash -l
#
# this is the global fronts image transfer task
#

THE_DATE=${1:-`date --date yesterday +%Y%m%d`}

LOCAL_DIR="/scratch2/NCEPDEV/stmp1/Samira.Ardani/images/class-4/fronts"
REMOTE_DIR="/home/www/polar/global/fronts/archive/images"

REMOTE_ID='emc.rtofs@emcrzdm.ncep.noaa.gov'

echo "Starting global fronts image transfer at `date`"

ssh ${REMOTE_ID} mkdir -p $REMOTE_DIR/$THE_DATE
ssh ${REMOTE_ID} mkdir -p $REMOTE_DIR/js
scp -r -q -o "BatchMode yes" $LOCAL_DIR/${THE_DATE} ${REMOTE_ID}:$REMOTE_DIR/.
scp -r -q -o "BatchMode yes" $LOCAL_DIR/stats ${REMOTE_ID}:$REMOTE_DIR/.
scp -r -q -o "BatchMode yes" $LOCAL_DIR/blockDates.js ${REMOTE_ID}:$REMOTE_DIR/.

echo "Completed global fronts image transfer at `date`"

# cleanup
rm -rf $LOCAL_DIR/${THE_DATE}

exit
