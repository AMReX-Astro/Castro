#!/bin/sh

# loop through all of the directories and do a make clean, to remove
# build temporaries

top=`pwd`

for p in `find . -mindepth 2 -maxdepth 2 -type d`
  do
    cd $p
    make realclean > /dev/null
    cd ${top}
  done

