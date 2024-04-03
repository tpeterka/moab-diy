#include "main.hpp"
#include "opts.h"

int main(int argc, char**argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    // command line options
    int     dim = 3;                            // domain dimensionality
    bool    help;

    using namespace opts;
    Options ops;
    ops
        >> Option('d', "dimension", dim,            "domain dimensionality")
        >> Option('h', "help",      help,           "show help")
        ;

    if (!ops.parse(argc,argv) || help)
    {
        if (world.rank() == 0)
        {
            std::cout << "Usage: " << argv[0] << " [OPTIONS]\n";
            std::cout << "Imports a MOAB parallel partition into a DIY block decomposition\n";
            std::cout << ops;
        }
        return 1;
    }

    // moab options
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

#if 0

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

    NeighBlocks send_neigh_blocks, recv_neigh_blocks;   // map of (neighboring rank, vector of my neighboring blocks) to send, receive

    // loop over local blocks, filling in the map of entities to parts
    Range::iterator part_it;
    std::map<long, int> entity_part_map;                // entity global id -> part id
    for (part_it = parts.begin(); part_it != parts.end(); ++part_it )
    {
        int part_id;
        rval = mbi->tag_get_data(part_tag, &(*part_it), 1, &part_id); ERR;

        // get all entities in the block
        Range ents;
        rval = mbi->get_entities_by_handle(*part_it, ents); ERR;
        fmt::print(stderr, "part_id (block gid) = {} ents.size() = {}\n", part_id, ents.size());
        for (auto ents_it = ents.begin(); ents_it != ents.end(); ++ents_it)
        {
            int ent_global_id;                          // globally unique ID for the entity
            rval = mbi->tag_get_data(mbi->globalId_tag(), &(*ents_it), 1, &ent_global_id); ERR;
//             fmt::print(stderr, "Entity local id {} global id {}\n", mbi->id_from_handle(*ents_it), ent_global_id);
            entity_part_map.emplace(ent_global_id, part_id);
        }
    }

    // debug: print the entity_part_map
//     for (auto it = entity_part_map.begin(); it != entity_part_map.end(); it++)
//         fmt::print(stderr, "entity_part_map[{}] = {}\n", it->first, it->second);

    // loop over local blocks, filling in their entities, adding blocks to the master, assigner, creating their links
    for (part_it = parts.begin(); part_it != parts.end(); ++part_it )
    {
        // create the block
        Block* b = new Block;
        b->eh = *part_it;
        rval = mbi->tag_get_data(part_tag, &b->eh, 1, &b->gid); ERR;
//         fmt::print(stderr, "gid = {}\n", b->gid);

        // create a link for the block
        diy::Link*      link = new diy::Link;           // link is this block's neighborhood
        diy::BlockID    neighbor;                       // one neighboring block

        // get all entities in the block
        Range ents;
        rval = mbi->get_entities_by_handle(b->eh, ents); ERR;
//         fmt::print(stderr, "ents.size() = {}\n", ents.size());

        // for all entities in local block
        for (auto ents_it = ents.begin(); ents_it != ents.end(); ++ents_it)
        {
            // get vertices comprising the element
            Range verts;
            rval = mbi->get_connectivity(&(*ents_it), 1, verts); ERR;
            for (auto verts_it = verts.begin(); verts_it != verts.end(); ++verts_it)
            {
                // get elements sharing this vertex
                Range adjs;
                rval = mbi->get_adjacencies(&(*verts_it), 1, dim, false, adjs, Interface::UNION); ERR;
//                 cout << CN::EntityTypeName(mbi->type_from_handle(*ents_it)) << " " << mbi->id_from_handle(*ents_it) << " Vertex " << mbi->id_from_handle(verts[i]) <<
//                     " is adjacent to " << endl;
//                 adjs.print();

                // iterate over adjacent elements
                for (auto adjs_iter = adjs.begin(); adjs_iter != adjs.end(); adjs_iter++)
                {
                    long elem_global_id;
                    int neigh_gid;
                    rval = mbi->tag_get_data(mbi->globalId_tag(), &(*adjs_iter), 1, &elem_global_id); ERR;
                    auto map_it = entity_part_map.find(elem_global_id);
                    if (map_it == entity_part_map.end())
                    {
                        fmt::print(stderr, "Error: element global id {} not found in entity_part_map\n", elem_global_id);
                        abort();
                    }
                    else
                        neigh_gid = map_it->second;                                     // neighboring block containing the adjacent element

                    // add the block to the link
                    if (link->find(neigh_gid) == -1 && neigh_gid != b->gid)
                    {
                        neighbor.gid    = neigh_gid;
                        neighbor.proc   = world.rank();
                        link->add_neighbor(neighbor);
                    }
                }
            }
        }   // for all entities in the local block

        // add the block to the master
        master.add(b->gid, b, link);

        // add the block to the assigner
        assigner.set_rank(world.rank(), b->gid);

    }   // for all local blocks

    // get entities owned by my process and shared with other processors
    Range shared_ents;
    rval = pc->get_shared_entities(-1, shared_ents, MBVERTEX); ERR;

    // debug
//     cout << "shared entities: " << endl;
//     shared_ents.print();

    // shared entities owned by my process
    Range my_shared_ents;
    rval = pc->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &my_shared_ents ); ERR;

    // shared entities owned by other processes
    Range other_shared_ents;
    rval = pc->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_AND, -1, &other_shared_ents); ERR;

    // collect sending messages

    // for all vertices owned by my rank and shared by other ranks
    for (auto my_shared_it = my_shared_ents.begin(); my_shared_it != my_shared_ents.end(); my_shared_it++)
    {
        std::set<int> send_procs;                   // processes to which I will send
        pc->get_sharing_data(&(*my_shared_it), 1, send_procs);

        // debug
//         cout << CN::EntityTypeName(mbi->type_from_handle(*my_shared_it)) << " " << mbi->id_from_handle(*my_shared_it) <<
//             " will be sent to " << send_procs.size() <<  " procs:" << endl;
//         for (auto send_procs_iter = send_procs.begin(); send_procs_iter != send_procs.end(); send_procs_iter++)
//             fmt::print(stderr, "*send_procs_iter = {}\n", *send_procs_iter);

        // assemble a set of block gids for each unique process found

        // get elements sharing this vertex
        Range adjs;
        rval = mbi->get_adjacencies(&(*my_shared_it), 1, dim, false, adjs, Interface::UNION); ERR;

        // iterate over elements sharing the vertex
        for (auto adjs_iter = adjs.begin(); adjs_iter != adjs.end(); adjs_iter++)
        {
            // get block gid containing the element
            long elem_global_id;
            int neigh_gid;
            rval = mbi->tag_get_data(mbi->globalId_tag(), &(*adjs_iter), 1, &elem_global_id); ERR;
            auto map_it = entity_part_map.find(elem_global_id);
            if (map_it == entity_part_map.end())
            {
                fmt::print(stderr, "Error: element global id {} not found in entity_part_map\n", elem_global_id);
                abort();
            }
            else
                neigh_gid = map_it->second;                                     // neighboring block containing the adjacent element

            // for all procs sharing this vertex, add block gid to the map of messages to send
            for (auto send_procs_iter = send_procs.begin(); send_procs_iter != send_procs.end(); send_procs_iter++)
            {
                auto found_neigh = send_neigh_blocks.find(*send_procs_iter);
                if (found_neigh == send_neigh_blocks.end()) // rank does not exist yet
                {
                    std::pair<int, set<int>> new_neigh;
                    new_neigh.first = *send_procs_iter;
                    new_neigh.second.insert(neigh_gid);
                    send_neigh_blocks.insert(new_neigh);

                    // debug
//                     fmt::print(stderr, "inserting new new_neigh first {} second [{}]\n", new_neigh.first, fmt::join(new_neigh.second, ","));
                }
                else        // rank exists already
                {
                    found_neigh->second.insert(neigh_gid);

                    // debug
//                     fmt::print(stderr, "inserting found new_neigh first {} second [{}]\n", found_neigh->first, fmt::join(found_neigh->second, ","));
                }
            }   // procs sharing the vertex
        }   // elements sharing the vertex
    }   // for all vertices owned by my rank and shared by other ranks

    // debug: print send_neigh_blocks
    for (auto neigh_it = send_neigh_blocks.begin(); neigh_it != send_neigh_blocks.end(); neigh_it++)
        fmt::print(stderr, "sending neighboring blocks to proc {}: [{}]\n", neigh_it->first, fmt::join(neigh_it->second, ","));

    // TODO: post messages

    // collect receiving messages

    // for all shared vertices owned by other ranks
    for (auto other_shared_it = other_shared_ents.begin(); other_shared_it != other_shared_ents.end(); other_shared_it++)
    {
        std::set<int> recv_procs;                   // processes from which I will receive
        pc->get_sharing_data(&(*other_shared_it), 1, recv_procs);

        // debug
//         cout << CN::EntityTypeName(mbi->type_from_handle(*other_shared_it)) << " " << mbi->id_from_handle(*other_shared_it) <<
//             " will be received from " << recv_procs.size() <<  " procs:" << endl;
//         for (auto recv_procs_iter = recv_procs.begin(); recv_procs_iter != recv_procs.end(); recv_procs_iter++)
//             fmt::print(stderr, "{}\n", *recv_procs_iter);

        // assemble a set of block gids for each unique process found

        // get elements sharing this vertex
        Range adjs;
        rval = mbi->get_adjacencies(&(*other_shared_it), 1, dim, false, adjs, Interface::UNION); ERR;

        // iterate over elements sharing the vertex
        for (auto adjs_iter = adjs.begin(); adjs_iter != adjs.end(); adjs_iter++)
        {
            // get block gid containing the element
            long elem_global_id;
            int neigh_gid;
            rval = mbi->tag_get_data(mbi->globalId_tag(), &(*adjs_iter), 1, &elem_global_id); ERR;
            auto map_it = entity_part_map.find(elem_global_id);
            if (map_it == entity_part_map.end())
            {
                fmt::print(stderr, "Error: element global id {} not found in entity_part_map\n", elem_global_id);
                abort();
            }
            else
                neigh_gid = map_it->second;                                     // neighboring block containing the adjacent element

            // for all procs sharing this vertex, add block gid to the map of messages to receive
            for (auto recv_procs_iter = recv_procs.begin(); recv_procs_iter != recv_procs.end(); recv_procs_iter++)
            {
                auto found_neigh = recv_neigh_blocks.find(*recv_procs_iter);
                if (found_neigh == recv_neigh_blocks.end()) // rank does not exist yet
                {
                    std::pair<int, set<int>> new_neigh;
                    new_neigh.first = *recv_procs_iter;
                    new_neigh.second.insert(neigh_gid);
                    recv_neigh_blocks.insert(new_neigh);

                    // debug
//                     fmt::print(stderr, "inserting new new_neigh first {} second [{}]\n", new_neigh.first, fmt::join(new_neigh.second, ","));
                }
                else        // rank exists already
                {
                    found_neigh->second.insert(neigh_gid);

                    // debug
//                     fmt::print(stderr, "inserting found new_neigh first {} second [{}]\n", found_neigh->first, fmt::join(found_neigh->second, ","));
                }
            }   // procs sharing the vertex
        }   // elements sharing the vertex
    }   // for all shared vertices owned by other ranks

    // debug: print recv_neigh_blocks
    for (auto neigh_it = recv_neigh_blocks.begin(); neigh_it != recv_neigh_blocks.end(); neigh_it++)
        fmt::print(stderr, "receiving neighboring blocks from proc {}: [{}]\n", neigh_it->first, fmt::join(neigh_it->second, ","));


    // TODO post recv messages

    // TODO: exchange link info with remote blocks

    // TODO: update the link

    // debug: print the link
    //         fmt::print(stderr, "Link for block gid {} has size {}:\n", b->gid, link->size());
    //         for (auto i = 0; i < link->size(); i++)
    //             fmt::print(stderr, "[gid, proc] = [{}, {}]\n", link->target(i).gid, link->target(i).proc);

    // write output file for debugging
    rval = mbi->write_file(outfile.c_str(), 0, write_opts.c_str(), &root, 1); ERR;

    return 0;
}
