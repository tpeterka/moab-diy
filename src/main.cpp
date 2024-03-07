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
    std::string infile      = "/home/tpeterka/software/spack/var/spack/environments/moab-example-env/moab-example/sample_data/mpas_2d_source_p128.h5m";
    std::string read_opts   = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;DEBUG_IO=3;";
    std::string outfile     = "example1.h5m";
    std::string write_opts  = "PARALLEL=WRITE_PART;DEBUG_IO=6";

    // debug
    fmt::print(stderr, "producer: local comm rank {} size {} metadata {} passthru {}\n",
            local_.rank(), local_.size(), metadata, passthru);

    int nafc = 0;               // number of times after file close callback was called

    if (shared)                 // single process, MetadataVOL test
        fmt::print(stderr, "producer: using shared mode MetadataVOL plugin created by prod-con\n");
    else                        // normal multiprocess, DistMetadataVOL plugin
    {
        l5::DistMetadataVOL& vol_plugin = l5::DistMetadataVOL::create_DistMetadataVOL(local, intercomms);

        // set lowfive properties
        if (passthru)
        {
            // debug
            fmt::print(stderr, "*** producer setting passthru mode\n");

            vol_plugin.set_passthru(outfile, "*");
        }
        if (metadata)
        {
            // debug
            fmt::print(stderr, "*** producer setting memory mode\n");

            vol_plugin.set_memory(outfile, "*");
        }
        vol_plugin.set_passthru(infile, "*");      // infile comes from disk
        vol_plugin.set_intercomm(infile, "*", 0);

        vol_plugin.set_keep(true);
        vol_plugin.serve_on_close = false;

        // set a callback to broadcast/receive files before a file open
        vol_plugin.set_before_file_open([&](const std::string& name)
        {
            if (name == outfile)
                vol_plugin.broadcast_files();
        });

        // set a callback to serve files after a file close
        vol_plugin.set_after_file_close([&](const std::string& name)
        {
            fmt::print(stderr, "##### producer afc 1: name = {}\n", name);
            if (name != outfile)
                return;

            fmt::print(stderr, "##### producer afc 2: name = {} nafc = {}\n", name, nafc);
            if (local_.rank() == 0)
            {
                if (nafc > 0)
                {
                    if (!passthru)
                    {
                        vol_plugin.serve_all();
                        vol_plugin.serve_all();
                    }
                }
            }
            else
            {
                if (!passthru)
                {
                    vol_plugin.serve_all();
                    vol_plugin.serve_all();
                }
            }

            nafc++;
            fmt::print(stderr, "##### producer afc 3: name = {} nafc = {}\n", name, nafc);
        });
    }

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

#if 0

    // create mesh in memory
    fmt::print(stderr, "*** producer generating synthetic mesh in memory ***\n");
    PrepMesh(mesh_type, mesh_size, mesh_slab, mbi, pc, root, factor, false);
    fmt::print(stderr, "*** producer after creating mesh in memory ***\n");

#else

    // or

    // read file
    fmt::print(stderr, "*** producer reading input file ***\n");
    rval = mbi->load_file(infile.c_str(), &root, read_opts.c_str() ); ERR(rval);
    fmt::print(stderr, "*** producer after reading file ***\n");

#endif

    // write file
    rval = mbi->write_file(outfile.c_str(), 0, write_opts.c_str(), &root, 1); ERR(rval);

    // debug
    fmt::print(stderr, "*** producer after writing file ***\n");

    // signal the consumer that data are ready
    if (passthru && !metadata && !shared)
    {
        for (auto& intercomm: intercomms)
            diy_comm(intercomm).barrier();
    }
    fmt::print(stderr, "*** producer after barrier completed! ***\n");

    return 0;
}
