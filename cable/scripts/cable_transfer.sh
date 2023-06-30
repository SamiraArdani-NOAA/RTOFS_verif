#!/bin/bash -l
#
# this is the cable image transfer task
#

LOCAL_DIR="/scratch1/NCEPDEV/stmp2/Samira.Ardani/images/class-4/cable"
REMOTE_DIR="/home/www/polar/global/class-4/cable/archive/images"

REMOTE_ID='emc.rtofs@emcrzdm.ncep.noaa.gov'

echo "Starting cable image transfer at `date`"

cd $LOCAL_DIR

for THE_DATE in *; do
  ssh ${REMOTE_ID} mkdir -p $REMOTE_DIR/$THE_DATE/
  scp -r -q -o "BatchMode yes" $LOCAL_DIR/${THE_DATE} ${REMOTE_ID}:$REMOTE_DIR/.
  rm -rf $LOCAL_DIR/${THE_DATE}
done

echo "Completed cable image transfer at `date`"

exit
