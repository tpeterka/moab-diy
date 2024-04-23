# moab-diy
Importing MOAB decomposition into DIY

# Instructions for Building and Running

Installation is done through Spack. If you don't have Spack installed or if Spack is new to you, go [here](https://spack.readthedocs.io/en/latest/) first.

## Setting up Spack environment

### First time: create and load the Spack environment

```
git clone https://github.com/tpeterka/moab-diy
cd /path/to/moab-diy
source ./create-env.sh     # requires being in the same directory to work properly
```

### Subsequent times: load the Spack environment

```
source /path/to/moab-diy/load-env.sh
```

----

## Building moab-diy

```
cd build
rm CMakeCache.txt
cmake .. \
-DCMAKE_INSTALL_PREFIX=/path/to/moab-diy/install
make -j install
```

-----

## Running moab-diy

There are three options for input data. You can generate a synthetic mesh, either tetrahedral or hexahedral, or read in input file as follows.

3d tetrahedral synthetic mesh

```
cd build/src
mpiexec -n <procs> ./main -i tet -s <size per side, optional, 10 is default>
```

3d hexahedral synthetic mesh

```
cd build/src
mpiexec -n <procs> ./main -i hex -s <size per side, optional, 10 is default>
```

2d or 3d input file that has been pre-partitioned. Sample input files are located in the sample_data directory.

```
cd build/src
mpiexec -n <procs> ./main -i file -f filename
```

-----
