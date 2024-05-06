#!/bin/bash

export SPACKENV=moab-diy-env
export YAML=$PWD/env.yaml

# create spack environment
echo "creating spack environment $SPACKENV"
spack env deactivate > /dev/null 2>&1
spack env remove -y $SPACKENV > /dev/null 2>&1
spack env create $SPACKENV $YAML

# activate environment
echo "activating spack environment"
spack env activate $SPACKENV

# install everything in environment
echo "installing dependencies in environment"
spack install

# reset the environment (workaround for spack behavior)
spack env deactivate
spack env activate $SPACKENV

# set build flags
echo "setting flags for building moab-example"
export DIY_PATH=`spack location -i diy`
export FMT_PATH=`spack location -i fmt`
export MOAB_PATH=`spack location -i moab`
export HDF5_PATH=`spack location -i hdf5`

# set LD_LIBRARY_PATH
echo "setting flags for running moab-example"
export LD_LIBRARY_PATH=$MOAB_PATH/lib:$LD_LIBRARY_PATH

