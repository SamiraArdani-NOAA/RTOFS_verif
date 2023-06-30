#!/bin/ksh

module purge
module use -a /contrib/anaconda/modulefiles
module load intel/2022.1.2 impi/2022.1.2
module load ncl/6.6.2
module load imagemagick/7.0.8-53
module load anaconda/latest
module list

echo "Started at `date`"

# cut date is one day before the run date
RunDate=${1:-$(date +%Y%m%d)}
CutDate=$(date --date="${RunDate} -1 day" +%Y%m%d)

font=${2:-10}

baseDir='/scratch1/NCEPDEV/marine/Zulema.Garraffo/rtofs_profile_plots'
tmpDir="/scratch1/NCEPDEV/stmp2/Samira.Ardani/RTOFS-DA"
dataDir="${baseDir}/rtofs.${RunDate}/ncoda/logs/profile_qc"
NEWHOME="/scratch1/NCEPDEV/stmp2/Samira.Ardani"
srcDir="${NEWHOME}/github/RTOFS_verif/ncoda/frames"
polarDir="/home/www/polar/develop/OpenLayers/profiles/${CutDate}"

# create tempdir
mkdir -p $tmpDir/png  $tmpDir/ps

# extract images from the four gmeta files of interest
for platform in xbt argo glider animal tesac buoy; do
  echo "Converting $platform gmeta to postscript"
  ctrans -device ps.color -simulatebg -font $font > ${tmpDir}/ps/gmeta_${platform}.tmp ${dataDir}/prof_${platform}_qc.${CutDate}00.gmeta
  cd ${tmpDir}/ps
  psplit gmeta_${platform}.tmp $platform
  cd ..
  
done  

# convert ps to png
echo 'Converting postscript to png'
cd ${tmpDir}/ps
maxtasks=20
ntasks=0
for file in *.ps; do
  ntasks=$(($ntasks + 1))
  file=$(basename $file '.ps')
  #echo $file
  convert -density 400 ${file}.ps -resize 800 -type Palette ${tmpDir}/png/${file}.png &
  if (( $ntasks >= $maxtasks )); then
    wait
    ntasks=0
  fi
done
wait
cd $tmpDir

# parse the output file
echo 'Parsing the out file to create markers.js'
python ${srcDir}/ush/profile_parser.py ${dataDir}/prof_qc.${CutDate}00.out
rm -f ${tmpDir}/png/buoy*.png

# upload to polar
echo 'Uploading to polar'
EMCRZDM='140.90.100.206'
USER='emc.rtofs'
ssh ${USER}@${EMCRZDM} mkdir -p ${polarDir}
scp ${tmpDir}/markers.js ${USER}@${EMCRZDM}:${polarDir}/.
scp -q ${tmpDir}/png/* ${USER}@${EMCRZDM}:${polarDir}/.

# cleanup
#rm -rf ${tmpDir}/*

echo "Completed at `date`" 
