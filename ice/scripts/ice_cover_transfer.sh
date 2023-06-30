#!/bin/bash -l
#
# this is the ice image transfer task
#

THE_DATE=${1:-`date --date yesterday +%Y%m%d`}

LOCAL_DIR="/scratch1/NCEPDEV/stmp2/Samira.Ardani/images/ice"
REMOTE_DIR="/home/www/polar/global/ice/archive/images"

REMOTE_ID='emc.rtofs@emcrzdm.ncep.noaa.gov'

echo "Starting ice image transfer at `date`"

ssh ${REMOTE_ID} mkdir -p $REMOTE_DIR/$THE_DATE
scp -r -q -o "BatchMode yes" $LOCAL_DIR/${THE_DATE} ${REMOTE_ID}:$REMOTE_DIR/.

echo "Completed ice image transfer at `date`"

# cleanup
rm -rf $LOCAL_DIR/${THE_DATE}

exit
