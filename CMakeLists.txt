project                     (moab-diy)
cmake_minimum_required      (VERSION 3.12)

option                      (wrapped_mpi    "MPI compiler wrapper requires no further MPI libs"     OFF)

# C++11
if                          (NOT DEFINED CMAKE_CXX_STANDARD)
    set                     (CMAKE_CXX_STANDARD 11)
endif                       ()

set                         (CMAKE_POSITION_INDEPENDENT_CODE ON)

# ignore variadic args warning (no other good solution)
set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-varargs")

# Default to Release
if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property            (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)

# MPI
if                          (NOT wrapped_mpi)
    find_package            (MPI REQUIRED)
    set                     (libraries MPI::MPI_C)
endif                       ()

# Threads
find_package                (Threads)

# DIY
if                          (NOT DEFINED DIY_PATH AND DEFINED ENV{DIY_PATH})
    set                     (DIY_PATH $ENV{DIY_PATH})
endif                       ()

if                          (DIY_PATH)
    message                 ("Found DIY_PATH")
    set                     (DIY_INCLUDE_DIR ${DIY_PATH}/include)
else                        ()
    message                 ("Looking for DIY")
    find_package            (diy REQUIRED)
endif                       ()
message                     ("DIY_INCLUDE_DIR = ${DIY_INCLUDE_DIR}")
include_directories         (SYSTEM ${DIY_INCLUDE_DIR})

# fmt
find_package                (fmt REQUIRED)

# Moab
if                          (NOT DEFINED MOAB_PATH AND DEFINED ENV{MOAB_PATH})
    set                     (MOAB_PATH $ENV{MOAB_PATH})
endif                       ()

if                          (MOAB_PATH)
    message                 ("Found MOAB_PATH")
    set                     (MOAB_LIBRARY ${MOAB_PATH}/lib/libMOAB.so)
    set                     (MOAB_IMESH_LIBRARY ${MOAB_PATH}/lib/libiMesh.so)
    set                     (MOAB_INCLUDE_DIR ${MOAB_PATH}/include)
    add_definitions         (-DMOAB_PATH)
else                        ()
    message                 ("Looking for Moab")
    find_package            (moab REQUIRED)
endif                       ()
message                     ("MOAB_LIBRARY = ${MOAB_LIBRARY}")
message                     ("MOAB_IMESH_LIBRARY = ${MOAB_IMESH_LIBRARY}")
message                     ("MOAB_INCLUDE_DIR = ${MOAB_INCLUDE_DIR}")
set                         (libraries ${libraries} ${MOAB_LIBRARY} ${MOAB_IMESH_LIBRARY})
include_directories         (SYSTEM ${MOAB_INCLUDE_DIR})

# HDF5
if                          (NOT DEFINED HDF5_PATH AND DEFINED ENV{HDF5_PATH})
    set                     (HDF5_PATH $ENV{HDF5_PATH})
endif                       ()

if                          (HDF5_PATH)
    set                     (HDF5_LIBRARY ${HDF5_PATH}/lib/libhdf5.so)
    set                     (HDF5_HL_LIBRARY ${HDF5_PATH}/lib/libhdf5_hl.so)
    set                     (HDF5_HL_INCLUDE_DIR ${HDF5_PATH}/include)
  message                   ("HDF5_LIBRARY = " ${HDF5_LIBRARY})
  message                   ("HDF5_HL_LIBRARY = " ${HDF5_HL_LIBRARY})
  message                   ("HDF5_INCLUDE_DIR = " ${HDF5_INCLUDE_DIR})
  message                   ("Z_LIBRARY = " ${Z_LIBRARY})
  set                       (HDF5_LIBRARY ${HDF5_HL_LIBRARY} ${HDF5_LIBRARY} ${Z_LIBRARY})
else                        ()
  message                   ("Looking for HDF5")
  find_package              (HDF5 REQUIRED COMPONENTS C HL)
  message                   ("HDF5_C_LIBRARIES = " ${HDF5_C_LIBRARIES})
  message                   ("HDF5_HL_LIBRARIES = " ${HDF5_HL_LIBRARIES})
  message                   ("HDF5_C_LIBRARY_hdf5 = " ${HDF5_C_LIBRARY_hdf5})
  message                   ("HDF5_C_LIBRARY_z = " ${HDF5_C_LIBRARY_z})
  message                   ("HDF5_INCLUDE_DIR = " ${HDF5_INCLUDE_DIR})
  set                       (HDF5_LIBRARY ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARY_hdf5} ${HDF5_C_LIBRARY_z})
endif                       ()
set                         (libraries ${libraries} ${HDF5_LIBRARY})
include_directories         (SYSTEM ${HDF5_INCLUDE_DIR})

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${CMAKE_CURRENT_SOURCE_DIR}/include
                            SYSTEM ${MPI_C_INCLUDE_PATH})

# Libraries
set                         (libraries
                             ${libraries}
                             ${CMAKE_THREAD_LIBS_INIT}
                             ${CMAKE_DL_LIBS})

add_subdirectory            (src)

