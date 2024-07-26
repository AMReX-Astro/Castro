#!/bin/sh

# clone Castro and all associated repos, then output the necessary
# environment variables to build with this source.

branch=master

amrex_server=https://github.com/AMReX-Codes
astro_server=https://github.com/AMReX-Astro
micro_server=https://github.com/StarKiller-astro


#-----------------------------------------------------------------------------

pwd=`pwd`

# Castro
echo "cloning Castro"
git clone ${astro_server}/Castro.git
cd Castro
git checkout ${branch}
cd ..


# AMReX
echo " "
echo "cloning amrex"
git clone ${amrex_server}/amrex.git
cd amrex
git checkout ${branch}
cd ..


# Microphysics
echo " "
echo "cloning Microphysics"
git clone ${micro_server}/Microphysics.git
cd Microphysics
git checkout ${branch}
cd ..


# output the necessary environment changes
if [ -f castro_exports.sh ]; then
    echo "removing old castro_exports.sh"
    rm -f castro_exports.sh
fi

cat >> castro_exports.sh << EOF
export CASTRO_HOME="${pwd}/Castro"
export MICROPHYSICS_HOME="${pwd}/Microphysics"
export AMREX_HOME="${pwd}/amrex"
EOF

echo " "
echo "source castro_exports.sh to setup the environment for building"




