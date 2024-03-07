#include "main.hpp"

herr_t fail_on_hdf5_error(hid_t stack_id, void*)
{
    H5Eprint(stack_id, stderr);
    fmt::print(stderr, "An HDF5 error was detected. Terminating.\n");
    exit(1);
}

int main(int argc, char**argv)
{
    diy::mpi::communicator local_(local);
    std::string infile      = "/home/tpeterka/software/moab-diy/sample_data/mpas_2d_source_p128.h5m";
    std::string read_opts   = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;DEBUG_IO=3;";

    // create moab mesh
    int                             mesh_type = 0;                          // source mesh type (0 = hex, 1 = tet)
    int                             mesh_size = 10;                        // source mesh size per side
    int                             mesh_slab = 0;                          // block shape (0 = cubes; 1 = slabs)
    double                          factor = 1.0;                           // scaling factor on field values
    Interface*                      mbi = new Core();                       // moab interface
    ParallelComm*                   pc  = new ParallelComm(mbi, local);     // moab communicator
    EntityHandle                    root;
    ErrorCode                       rval;
    rval = mbi->create_meshset(MESHSET_SET, root); ERR(rval);

#if 1

    // create mesh in memory
    PrepMesh(mesh_type, mesh_size, mesh_slab, mbi, pc, root, factor, false);

#else

    // or

    // read file
    rval = mbi->load_file(infile.c_str(), &root, read_opts.c_str() ); ERR(rval);

#endif

    return 0;
}
