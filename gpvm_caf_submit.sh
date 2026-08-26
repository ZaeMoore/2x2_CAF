#!/bin/bash  

USERNAME=${zaemoore}

TOPDIR="/pnfs/dune/scratch/users/zaemoore"

/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer shell --shell=/bin/bash \
-B /cvmfs,/exp,/nashome,/pnfs/dune,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf --ipc --pid \
/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

#Taken from ND_CAFMaker, some of these might not be necessary
setup cmake v3_22_2
setup gcc v9_3_0
setup pycurl
setup ifdhc
setup geant4 v4_11_0_p01c -q e20:debug
setup dk2nugenie   v01_10_01k -q debug:e20
setup genie_xsec   v3_04_00 -q AR2320i00000:e1000:k250
setup genie_phyopt v3_04_00 -q dkcharmtau
setup jobsub_client
setup eigen v3_3_5
setup duneanaobj v03_01_00 -q e20:prof
setup srproxy v00.43 -q py3913
setup hdf5 v1_10_5a -q e20
setup fhiclcpp v4_15_03 -q debug:e20

# edep-sim needs to know where a certain GEANT .cmake file is...
G4_cmake_file=`find ${GEANT4_FQ_DIR}/lib -name 'Geant4Config.cmake'`
export Geant4_DIR=`dirname $G4_cmake_file`

# edep-sim needs to have the GEANT bin directory in the path
export PATH=$PATH:$GEANT4_FQ_DIR/bin

# shut up ROOT include errors
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$GENIE_INC/GENIE

# Add pyGeoEff to pythonpath
export PYTHONPATH=${PYTHONPATH}:${PWD}/DUNE_ND_GeoEff/lib/

# duneananobj needs to be in the libs too
export LD_LIBRARY_PATH=${DUNEANAOBJ_LIB}:$LD_LIBRARY_PATH

# finally, add our lib & bin to the paths
mydir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export LD_LIBRARY_PATH=$mydir/lib:$LD_LIBRARY_PATH
export PATH=$mydir/bin:$PATH

#cp minirun5_flat.list

