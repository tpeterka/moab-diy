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

In src/main.cpp, you can control whether to create a synthetic mesh in memory or read an input file with
`#if 1` or `#if 0` at/around line 384.

```
cd build/src
mpiexec -n <procs> ./main -d <dimension of elements, usually 2 or 3, defaults to 3>
```

Note that the sample data file is 2d, and requires `-d 2` in the command line

-----
