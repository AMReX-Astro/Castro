#!/bin/sh -f

if [ ! "$1" ]; then
  echo "usage: chainsub.sh jobid number"
  exit -1
fi

if [ ! "$2" ]; then
  echo "usage: chainsub.sh jobid number"
  exit -1
fi

echo chaining $2 jobs starting with $1

oldjob=$1

for count in `seq 1 1 $2`
do
  echo starting job $count to depend on $oldjob
  aout=`qsub -W depend=afterany:${oldjob} titan.run`
  echo "   " jobid: $aout
  echo " "
  oldjob=$aout
done
