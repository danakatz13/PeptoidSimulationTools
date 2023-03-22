#!/bin/bash -f

if [ -f groupfile ];
then
  rm groupfile
fi

nrep=`wc temperatures.dat | awk '{print $1}'`
echo $nrep
count=0
for TEMP in `cat temperatures.dat`
do
  let COUNT+=1
  REP=`printf "%03d" $COUNT`
  echo "TEMPERATURE: $TEMP K ==> FILE: prod.in.$REP"
  sed "s/XXXXX/$TEMP/g" prod.in > temp
  sed "s/RANDOM_NUMBER/$RANDOM/g" temp > remd.mdin.$REP
  echo "-O -rem 1 -remlog rem.log -i remd.mdin.$REP -o prod.out.$REP -c equilibrate.rst.$REP -r remd.rst.$REP -x remd.mdcrd.$REP -inf remd.mdinfo.$REP -p design_1.parm7" >> remd.groupfile

  rm -f temp
done
echo "#" >> groupfile

echo "N REPLICAS  = $nrep"
echo " Done."
