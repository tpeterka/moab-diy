# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class MoabDiy(CMakePackage):
    """Importing MOAB Decompostion into DIY."""

    homepage = "https://github.com/tpeterka/moab-diy.git"
    url      = "https://github.com/tpeterka/moab-diy.git"
    git      = "https://github.com/tpeterka/moab-diy.git"

    version('main', branch='main')

    depends_on('mpich')
    depends_on('hdf5+mpi+hl', type='link')
    depends_on('moab', type='link')
    depends_on('diy')
    depends_on('fmt')

    def cmake_args(self):
        args = ['-DCMAKE_C_COMPILER=%s' % self.spec['mpich'].mpicc,
                '-DCMAKE_CXX_COMPILER=%s' % self.spec['mpich'].mpicxx,
                '-DDIY_PATH=%s' % self.spec['diy'].prefix,
                '-DMOAB_PATH=%s' % self.spec['moab'].prefix,
                '-DFMT_PATH=%s' % self.spec['fmt'].prefix,
                '-DHDF5_PATH=%s' % self.spec['hdf5'].prefix,
                '-DBUILD_SHARED_LIBS=false']
        return args
