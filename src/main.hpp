#pragma once

#include    <vector>
#include    <cassert>

// diy
#include    <diy/master.hpp>
#include    <diy/decomposition.hpp>
#include    <diy/assigner.hpp>

// fmt
#include    <fmt/format.h>

// moab
#include    "iMesh.h"
#include    "MBiMesh.hpp"
#include    "moab/Core.hpp"
#include    "moab/Range.hpp"
#include    "MBTagConventions.hpp"
#include    "moab/ParallelComm.hpp"
#include    "moab/HomXform.hpp"
#include    "moab/ReadUtilIface.hpp"
#include    "moab/CN.hpp"

using mpi_comm      = MPI_Comm;
using diy_comm      = diy::mpi::communicator;
using Bounds        = diy::DiscreteBounds;

struct Block
{
  static void*    create()                                    { return new Block; }
  static void     destroy(void* b)                            { delete static_cast<Block*>(b); }
  static void     save(const void* b, diy::BinaryBuffer& bb)  { diy::save(bb, *static_cast<const Block*>(b)); }
  static void     load(void* b, diy::BinaryBuffer& bb)        { diy::load(bb, *static_cast<Block*>(b)); }

  // TODO: add user-defined member functions

  // TODO: add  user-defined data members
  int           gid;    // global id of this block
  EntityHandle  eh;     // root of everything moab-related in this block
};

#define ERR {if(rval!=MB_SUCCESS)printf("MOAB error at line %d in %s\n", __LINE__, __FILE__);}

using namespace moab;
using namespace std;

void hex_mesh_gen(int *mesh_size, Interface *mbint, EntityHandle *mesh_set,
        ParallelComm *mbpc, diy::RegularDecomposer<Bounds>& decomp, diy::RoundRobinAssigner& assign);
void tet_mesh_gen(int *mesh_size, Interface *mbint, EntityHandle *mesh_set,
        ParallelComm *mbpc, diy::RegularDecomposer<Bounds>& decomp, diy::RoundRobinAssigner& assign);
void create_hexes_and_verts(int *mesh_size, Interface *mbint, EntityHandle *mesh_set,
        diy::RegularDecomposer<Bounds>& decomp, diy::RoundRobinAssigner& assign, ParallelComm* mbpc);
void create_tets_and_verts(int *mesh_size, Interface *mbint, EntityHandle *mesh_set,
        diy::RegularDecomposer<Bounds>& decomp, diy::RoundRobinAssigner& assign, ParallelComm* mbpc);
void resolve_and_exchange(Interface *mbint, EntityHandle *mesh_set, ParallelComm *mbpc);

double PhysField(double x, double y, double z, double factor);

void PutElementField(Interface *mbi, EntityHandle eh, const char *tagname, double factor);

void GetElementField(Interface *mbi, EntityHandle eh, const char *tagname, double factor, MPI_Comm comm, bool debug);

void PutVertexField(Interface *mbi, EntityHandle eh, const char *tagname, double factor);

void GetVertexField(Interface *mbi, EntityHandle eh, const char *tagname, double factor, MPI_Comm comm, bool debug);

void PrintMeshStats(Interface *mbint, EntityHandle *mesh_set, ParallelComm *mbpc);

void PrepMesh(int src_type, int src_size, int slab, Interface* mbi, ParallelComm* pc, EntityHandle root, double factor,
        bool debug);
