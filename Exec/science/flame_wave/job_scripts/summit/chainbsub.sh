#!/bin/sh -f

if [ ! "$1" ]; then
  echo "usage: chainbsub.sh jobid number"
  exit -1
fi

if [ ! "$2" ]; then
  echo "usage: chainbsub.sh jobid number"
  exit -1
fi

echo chaining $2 jobs starting with $1

oldjob=$1

for count in `seq 1 1 $2`
do
  echo starting job $count to depend on $oldjob
  aout=`bsub -w "ended(${oldjob})" summit_gpu.submit`
  id=`echo $aout | head -n1 | cut -d'<' -f2 | cut -d'>' -f1`
  echo "   " jobid: $id
  echo " "
  oldjob=$id
  sleep 1
done
