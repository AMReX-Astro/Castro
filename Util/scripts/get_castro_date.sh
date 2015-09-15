#!/bin/sh

# clone Castro and all associated repos and reset them to a particular
# date/time.  Then output the necessary environment variables to build
# with this source.


server=https://github.com/BoxLib-Codes

#-----------------------------------------------------------------------------

date=$1

pwd=`pwd`

# MAESTRO
echo "cloning Castro"
git clone ${server}/Castro.git

echo " "
echo "resetting to before ${date}"
cd Castro
hash=`git rev-list -n 1 --before="$date" master`
git reset --hard ${hash}

cd ..

# BoxLib
echo " "
echo "cloning BoxLib"
git clone ${server}/BoxLib.git

echo " "
echo "resetting to before ${date}"
cd BoxLib
hash=`git rev-list -n 1 --before="$date" master`
git reset --hard ${hash}


# output the necessary environment changes
if [ -f exports.sh ]; then
    echo "removing old exports.sh"
    rm -f exports.sh
fi

cat >> exports.sh << EOF 
export CASTRO_DIR="${pwd}/Castro"
export BOXLIB_HOME="${pwd}/BoxLib"
EOF

echo " "
echo "source exports.sh to setup the environment for building"




