#include "main.hpp"

int main(int argc, char**argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    std::string infile      = "/home/tpeterka/software/moab-diy/sample_data/mpas_2d_source_p128.h5m";
    std::string read_opts   = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;DEBUG_IO=0;";
    std::string outfile     = "outfile.h5m";
    std::string write_opts  = "PARALLEL=WRITE_PART;DEBUG_IO=0";

    // create moab mesh
    int                             mesh_type = 1;                          // source mesh type (0 = hex, 1 = tet)
    int                             mesh_size = 10;                         // source mesh size per side
    int                             mesh_slab = 0;                          // block shape (0 = cubes; 1 = slabs)
    double                          factor = 1.0;                           // scaling factor on field values
    Interface*                      mbi = new Core();                       // moab interface
    ParallelComm*                   pc  = new ParallelComm(mbi, world);     // moab communicator
    EntityHandle                    root;
    ErrorCode                       rval;
    rval = mbi->create_meshset(MESHSET_SET, root); ERR;

#if 1

    // create mesh in memory
    PrepMesh(mesh_type, mesh_size, mesh_slab, mbi, pc, root, factor, false);

#else

    // or

    // read file
    rval = mbi->load_file(infile.c_str(), &root, read_opts.c_str() ); ERR;

#endif

    // query number of local parts in the parallel moab partition
    Range   parts;
    int     *part_values;
    Tag     part_tag;
    int     nblocks;            // local number of blocks = local number of parts
    int     tot_nblocks;        // total global number of blocks
    rval = mbi->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, part_tag); ERR;
    rval = mbi->get_entities_by_type_and_tag(0, MBENTITYSET, &part_tag, NULL, 1, parts); ERR;
    nblocks = parts.size();

    // sum all local parts to get global total number of blocks
    diy::mpi::all_reduce(world, nblocks, tot_nblocks, std::plus<int>());
    fmt::print(stderr, "nblocks = {} tot_nblocks = {}\n", nblocks, tot_nblocks);

    // master is in charge of local blocks
    diy::FileStorage          storage("./DIY.XXXXXX");
    diy::Master               master(world,             // MPI world communicator
                                     1,                 // # blocks to execute concurrently (1 means sequence through local blocks)
                                     -1,                // # blocks to keep in memory (-1 means keep all blocks in memory)
                                     &Block::create,    // block create function
                                     &Block::destroy,   // block destroy function
                                     &storage,          // storage location for out-of-core blocks (optional)
                                     &Block::save,      // block save function for out-of-core blocks (optional)
                                     &Block::load);     // block load function for out-of-core blocks (otpional)

    // create a dynamic assinger because we need to explicitly specify how many each rank has
    diy::DynamicAssigner assigner(world, world.size(), tot_nblocks);

    // loop over local blocks, filling in their entities, adding blocks to the master, assigner, creating their links
    Range::iterator part_it;
    for (part_it = parts.begin(); part_it != parts.end(); ++part_it )
    {
        // create the block
        Block* b = new Block;
        b->eh = *part_it;
        rval = mbi->tag_get_data(part_tag, &b->eh, 1, &b->gid);
        fmt::print(stderr, "gid = {}\n", b->gid);

        // create a link for the block
        diy::Link*      link = new diy::Link;           // link is this block's neighborhood
        diy::BlockID    neighbor;                       // one neighboring block
        // TODO fill in the link

        // get all entities in the block
        Range ents;
        rval = mbi->get_entities_by_handle(b->eh, ents); ERR;
        fmt::print(stderr, "ents.size() = {}\n", ents.size());
        for (Range::iterator ents_it = ents.begin(); ents_it != ents.end(); ++ents_it)
        {
            // get vertices comprising the element
            const EntityHandle* verts;
            int num_verts;
            rval = mbi->get_connectivity(*ents_it, verts, num_verts); ERR;
            for( int i = 0; i < num_verts; i++ )
            {
                // get block containing each vertex
                Range vert_sets;
                rval = mbi->get_adjacencies(&verts[i], 1, MBENTITYSET, false, vert_sets, Interface::UNION); ERR;
                cout << CN::EntityTypeName(mbi->type_from_handle(verts[i])) << " " << mbi->id_from_handle(verts[i]) << " is contained in sets: " << endl;
                vert_sets.print();

                // get elements sharing this vertex
                Range adjs;
                rval = mbi->get_adjacencies(&verts[i], 1, 3, false, adjs, Interface::UNION); ERR;
                cout << CN::EntityTypeName(mbi->type_from_handle(*ents_it)) << " " << mbi->id_from_handle(*ents_it) << " Vertex " << mbi->id_from_handle(verts[i]) <<
                    " is adjacencent to " << endl;
                adjs.print();

                // get block containing each element
                for (auto adj_it = adjs.begin(); adj_it != adjs.end(); ++adj_it)
                {
                    Range sets;
                    rval = mbi->get_adjacencies(&(*adj_it), 1, MBENTITYSET, false, sets, Interface::UNION); ERR;
                    cout << CN::EntityTypeName(mbi->type_from_handle(*adj_it)) << " " << mbi->id_from_handle(*adj_it) << " is contained in sets: " << endl;
                    sets.print();
                }
            }

//                 cout << CN::EntityTypeName(mbi->type_from_handle(*it)) << " " << mbi->id_from_handle(*it) << " vertex connectivity is: " << connect.print() << endl;
//                 fmt::print(stderr, "{} {} vertex connectivity is: ", CN::EntityTypeName(mbi->type_from_handle(*it)), mbi->id_from_handle(*it));
//                 std::vector<EntityHandle> parents;
//                 for (int i = 0; i < num_connect; i++)
//                 {
//                     parents.clear();
//                     rval = mbi->get_parent_meshsets(connect[i], parents); ERR;
//                     if (parents.size())
//                         fmt::print(stderr, "{} parents [{}] ", mbi->id_from_handle(connect[i]), fmt::join(parents, ","));
//                     fmt::print(stderr, "{} ", mbi->id_from_handle(connect[i]));
//                 }
//                 fmt::print(stderr, "\n");
//             }
        }

        // add the block to the master
        master.add(b->gid, b, link);

        // add the block to the assigner
        assigner.set_rank(world.rank(), b->gid);
    }

    // write output file for debugging
    rval = mbi->write_file(outfile.c_str(), 0, write_opts.c_str(), &root, 1); ERR;

    return 0;
}
