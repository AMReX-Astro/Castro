#!/bin/sh

# clone Castro and all associated repos and reset them to a particular
# date/time.  Then output the necessary environment variables to build
# with this source.


server=https://github.com/BoxLib-Codes
branch=development

#-----------------------------------------------------------------------------

date=$1

pwd=`pwd`

# Castro
echo "cloning Castro"
git clone ${server}/Castro.git

echo " "
echo "resetting to before ${date}"
cd Castro
git checkout ${branch}
hash=`git rev-list -n 1 --before="$date" $branch`
git reset --hard ${hash}

cd ..

# BoxLib
echo " "
echo "cloning BoxLib"
git clone ${server}/BoxLib.git

echo " "
echo "resetting to before ${date}"
cd BoxLib
git checkout ${branch}
hash=`git rev-list -n 1 --before="$date" $branch`
git reset --hard ${hash}

cd ..


# Microphysics
echo " "
echo "cloning Microphysics"
git checkout ${branch}
git clone ${server}/Microphysics.git

echo " "
echo "resetting to before ${date}"
cd Microphysics
hash=`git rev-list -n 1 --before="$date" $branch`
git reset --hard ${hash}

cd ..



# output the necessary environment changes
if [ -f exports.sh ]; then
    echo "removing old exports.sh"
    rm -f exports.sh
fi

cat >> exports.sh << EOF 
export CASTRO_HOME="${pwd}/Castro"
export MICROPHYSICS_HOME="${pwd}/Microphysics"
export BOXLIB_HOME="${pwd}/BoxLib"
EOF

echo " "
echo "source exports.sh to setup the environment for building"




