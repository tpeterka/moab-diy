#!/bin/bash

# activate the environment
export SPACKENV=moab-diy-env
spack env deactivate > /dev/null 2>&1
spack env activate $SPACKENV
echo "activated spack environment $SPACKENV"

echo "setting flags for building moab-example"
export MOAB_DIY_PATH=`spack location -i moab-example`
export DIY_PATH=`spack location -i diy`
export FMT_PATH=`spack location -i fmt`
export MOAB_PATH=`spack location -i moab`
export HDF5_PATH=`spack location -i hdf5`

echo "setting flags for running moab-example"
export LD_LIBRARY_PATH=$MOAB_PATH/lib:$LD_LIBRARY_PATH


