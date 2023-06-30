#!/bin/bash -l
#
# this is the satellite image transfer task
#

THE_DATE=${1:-`date --date="2 days ago" +%Y%m%d`}
#THE_DATE=${1:-`date --date yesterday +%Y%m%d`}

LOCAL_DIR="/scratch1/NCEPDEV/stmp2/Samira.Ardani/images/class-4/satellite"
REMOTE_DIR="/home/www/polar/global/class-1/archive/images"

REMOTE_ID='emc.rtofs@emcrzdm.ncep.noaa.gov'

echo "Starting Class-1 satellite for $THE_DATE image transfer at `date`"

ssh ${REMOTE_ID} mkdir -p $REMOTE_DIR/$THE_DATE
scp -r -q -o "BatchMode yes" $LOCAL_DIR/${THE_DATE} ${REMOTE_ID}:$REMOTE_DIR/.

echo "Completed Class-1 satellite for $THE_DATE image transfer at `date`"

# cleanup
rm -rf $LOCAL_DIR/${THE_DATE}

exit
