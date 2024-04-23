#include "main.hpp"
#include "opts.h"
#include "stdio.h"

// construct map of entities to parts (blocks)
void build_entity_part_map(const Range&     parts,                              // parts in the partition
                           Tag              part_tag,                           // moab part tag
                           Interface*       mbi,                                // moab interface
                           ParallelComm*    pc,                                 // moab parallel communicator
                           EntityPart&      entity_part_map)                    // (output) map of entity gid -> part (block) gid
{
    ErrorCode   rval;

    // loop over local blocks, filling in the map of entities to parts
    for (auto part_it = parts.begin(); part_it != parts.end(); ++part_it)
    {
        int part_id;
        rval = mbi->tag_get_data(part_tag, &(*part_it), 1, &part_id); ERR;

        // get all entities in the block
        Range ents;
        rval = mbi->get_entities_by_handle(*part_it, ents); ERR;

        // check first entity to see if global ids have been assigned; assign it if not
        int ent_global_id;
        rval = mbi->tag_get_data(mbi->globalId_tag(), &(*ents.begin()), 1, &ent_global_id); ERR;
        if (ent_global_id == -1)
        {
            // debug
            fmt::print(stderr, "No global ids were assigned; generating global ids now.\n");

            // determine 2d or 3d dimensionality from type of first element and assign global ids to elements
            EntityType type = mbi->type_from_handle(*(ents.begin()));
            if (type > MBVERTEX && type <= MBPOLYGON)
                pc->assign_global_ids(*part_it, 2);
            else if (type > MBPOLYGON && type <= MBPOLYHEDRON)
                pc->assign_global_ids(*part_it, 3);
            else
            {
                fmt::print(stderr, "Error: unknown entity type.\n");
                abort();
            }
            // also assign global ids to vertices
            pc->assign_global_ids(*part_it, 0);
        }

        // write the entity_part_map
        for (auto ents_it = ents.begin(); ents_it != ents.end(); ++ents_it)
        {
            rval = mbi->tag_get_data(mbi->globalId_tag(), &(*ents_it), 1, &ent_global_id); ERR;
            entity_part_map.emplace(ent_global_id, part_id);
        }
    }

    // debug: print the entity_part_map
//     for (auto it = entity_part_map.begin(); it != entity_part_map.end(); it++)
//         fmt::print(stderr, "entity_part_map[{}] = {}\n", it->first, it->second);
}

// loop over local blocks, filling in their entities, adding blocks to the master, assigner, creating their links
void collect_local_blocks(const Range&          parts,                          // parts in the moab partition
                          Tag                   part_tag,                       // moab part tag
                          Interface*            mbi,                            // moab interface
                          const EntityPart&     entity_part_map,                // map of entity gid -> part (block) gid
                          diy::Master&          master,                         // diy master
                          diy::DynamicAssigner& assigner)                       // diy assigner
{
    ErrorCode   rval;

    for (auto part_it = parts.begin(); part_it != parts.end(); ++part_it )
    {
        // create the block
        Block* b = new Block;
        b->eh = *part_it;
        rval = mbi->tag_get_data(part_tag, &b->eh, 1, &b->gid); ERR;

        // create a link for the block
        diy::Link*      link = new diy::Link;           // link is this block's neighborhood
        diy::BlockID    neighbor;                       // one neighboring block

        // get all entities in the block
        Range ents;
        rval = mbi->get_entities_by_handle(b->eh, ents); ERR;

        // for all entities in local block
        for (auto ents_it = ents.begin(); ents_it != ents.end(); ++ents_it)
        {
            // debug
//             int elem_global_id;
//             rval = mbi->tag_get_data(mbi->globalId_tag(), &(*ents_it), 1, &elem_global_id); ERR;
//             fmt::print(stderr, "elem_global_id = {}\n", elem_global_id);

            // get vertices comprising the element
            Range verts;
            rval = mbi->get_connectivity(&(*ents_it), 1, verts); ERR;
            for (auto verts_it = verts.begin(); verts_it != verts.end(); ++verts_it)
            {
                // debug
//                 int vert_global_id;
//                 rval = mbi->tag_get_data(mbi->globalId_tag(), &(*verts_it), 1, &vert_global_id); ERR;
//                 fmt::print(stderr, "vert_global_id = {}\n", vert_global_id);

                // get elements sharing this vertex
                // try 3d first, if empty, try 2d
                Range adjs;
                rval = mbi->get_adjacencies(&(*verts_it), 1, 3, false, adjs, Interface::UNION); ERR;
                if (adjs.size() == 0)
                    rval = mbi->get_adjacencies(&(*verts_it), 1, 2, false, adjs, Interface::UNION); ERR;
                if (adjs.size() == 0)
                {
                    fmt::print(stderr, "Error: no adjacencies found\n");
                    abort();
                }

                // iterate over adjacent elements
                for (auto adjs_it = adjs.begin(); adjs_it != adjs.end(); adjs_it++)
                {
                    int elem_global_id;
                    int neigh_gid;
                    rval = mbi->tag_get_data(mbi->globalId_tag(), &(*adjs_it), 1, &elem_global_id); ERR;
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
                        neighbor.proc   = master.communicator().rank();
                        link->add_neighbor(neighbor);
                    }
                }
            }
        }   // for all entities in the local block

        // add the block to the master
        master.add(b->gid, b, link);

        // add the block to the assigner
        assigner.set_rank(master.communicator().rank(), b->gid);

        // debug: print the link before remote exchange
//         fmt::print(stderr, "Link before exchange for block gid {} has size {}:\n", b->gid, link->size());
//         for (auto i = 0; i < link->size(); i++)
//             fmt::print(stderr, "[gid, proc] = [{}, {}]\n", link->target(i).gid, link->target(i).proc);

    }   // for all local blocks
}

// collect neighbor information for vertices shared across process boundaries
void collect_neigh_info(const Range&            shared_verts,                       // vertices shared by procs
                        Interface*              mbi,                                // moab interface
                        ParallelComm*           pc,                                 // moab parallel communicator
                        const EntityPart&       entity_part_map,                    // map of entity gid -> part (block) gid
                        ProcNeighBlocks&        shared_neigh_blocks)                // (output) neighbor info
{
    ErrorCode rval;

    // for all shared vertices
    for (auto shared_it = shared_verts.begin(); shared_it != shared_verts.end(); shared_it++)
    {
        // get set of procs with whom I share vertices
        std::set<int> shared_procs;                                               // processes with which I share entities
        pc->get_sharing_data(&(*shared_it), 1, shared_procs);

        // debug
//         fmt::print(stderr, "shared procs = [{}]\n", fmt::join(shared_procs, ","));
//         cout << CN::EntityTypeName(mbi->type_from_handle(*shared_it)) << " " << mbi->id_from_handle(*shared_it) <<
//             " will exchange with " << shared_procs.size() <<  " procs:" << endl;
//         for (auto shared_procs_iter = shared_procs.begin(); shared_procs_iter != shared_procs.end(); shared_procs_iter++)
//             fmt::print(stderr, "*shared_procs_iter = {}\n", *shared_procs_iter);

        // assemble a set of block gids for each unique process found

        // get elements sharing this vertex
        // try 3d first, if empty, try 2d
        Range adjs;
        rval = mbi->get_adjacencies(&(*shared_it), 1, 3, false, adjs, Interface::UNION); ERR;
        if (adjs.size() == 0)
            rval = mbi->get_adjacencies(&(*shared_it), 1, 2, false, adjs, Interface::UNION); ERR;
        if (adjs.size() == 0)
        {
            fmt::print(stderr, "Error: no adjacencies found\n");
            abort();
        }

        // iterate over elements sharing the vertex
        for (auto adjs_iter = adjs.begin(); adjs_iter != adjs.end(); adjs_iter++)
        {
            // get block gid containing the element
            int elem_global_id;
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

            int vert_global_id;                                                 // globally unique ID for the entity
            rval = mbi->tag_get_data(mbi->globalId_tag(), &(*shared_it), 1, &vert_global_id); ERR;

            // debug
//             cerr << " Vertex " << mbi->id_from_handle(*shared_it) << " global id " << vert_global_id << " neigh_gid " << neigh_gid << endl;

            // for all procs sharing this vertex, add block gid to the map of messages to send
            for (auto shared_procs_iter = shared_procs.begin(); shared_procs_iter != shared_procs.end(); shared_procs_iter++)
            {
                // moab sharing data includes my own proc, skip myself
                if (*shared_procs_iter == pc->rank())
                    continue;

                std::pair<int, set<int>> new_neigh;                             // (vertex gid, block gid)
                new_neigh.first = vert_global_id;
                new_neigh.second.insert(neigh_gid);

                auto found_proc = shared_neigh_blocks.find(*shared_procs_iter);
                if (found_proc == shared_neigh_blocks.end())                      // proc does not exist yet
                {
                    std::pair<int, NeighBlocks>  proc_new_neigh;                // (proc, new_neigh)
                    proc_new_neigh.first = *shared_procs_iter;
                    proc_new_neigh.second.insert(new_neigh);
                    shared_neigh_blocks.insert(proc_new_neigh);

                    // debug
//                     fmt::print(stderr, "shared_neigh_blocks: inserting new proc {} new new_neigh first {} second [{}]\n",
//                     proc_new_neigh.first, new_neigh.first, fmt::join(new_neigh.second, ","));
                }
                else                                                            // proc exists already
                {
                    auto found_neigh = found_proc->second.find(new_neigh.first);
                    if (found_neigh == found_proc->second.end())                // vertex gid does not exist yet
                    {
                        found_proc->second.insert(new_neigh);

                        // debug
//                         fmt::print(stderr, "shared_neigh_blocks: inserting existing proc {} new new_neigh first {} second [{}]\n",
//                                 found_proc.first, new_neigh.first, fmt::join(new_neigh.second, ","));
                    }
                    else                                                        // vertex gid exists already
                    {
                        found_neigh->second.insert(neigh_gid);

                        // debug
//                         fmt::print(stderr, "shared_neigh_blocks: inserting existing proc {} existing new_neigh first {} second [{}]\n",
//                                 found_proc.first, found_neigh.first, fmt::join(new_neigh.second, ","));
                    }
                }

            }   // procs sharing the vertex
        }   // elements sharing the vertex
    }   // for all vertices owned by my rank and shared by other ranks

    // debug: print shared_neigh_blocks
//     for (auto proc_it = shared_neigh_blocks.begin(); proc_it != shared_neigh_blocks.end(); proc_it++)
//     {
//         fmt::print(stderr, "sharing neighboring blocks with proc {}:\n", proc_it->first);
//         for (auto neigh_it = proc_it->second.begin(); neigh_it != proc_it->second.end(); neigh_it++)
//             fmt::print(stderr, "vid {}: block gids [{}]\n", neigh_it->first, fmt::join(neigh_it->second, ","));
//     }
}

// send neighbor info to other processes
void send_neigh_info(const ProcNeighBlocks&     shared_neigh_blocks,
                     diy::mpi::communicator&    comm)
{
    std::vector<int> msg;

    // for all procs in the neighbor info
    for (auto proc_it = shared_neigh_blocks.begin(); proc_it != shared_neigh_blocks.end(); proc_it++)
    {
        int dest_proc = proc_it->first;
        // for all vertex ids for this proc
        for (auto neigh_it = proc_it->second.begin(); neigh_it != proc_it->second.end(); neigh_it++)
        {
            msg.clear();
            msg.push_back(neigh_it->first);      // vertex gid

            // for all blocks for this vertex
            for (auto block_it = neigh_it->second.begin(); block_it != neigh_it->second.end(); block_it++)
                msg.push_back(*block_it);

            comm.send(dest_proc, 0, msg);

            // debug
//             fmt::print(stderr, "sent message to proc {}: [{}]\n", dest_proc, fmt::join(msg, ","));
        }
    }
}

// receive neighbor info from other processes
// update the links based on received info
void recv_neigh_info(const ProcNeighBlocks&     shared_neigh_blocks,        // neighbor blocks being shared
                     const ProcNeighBlocks&     extra_neigh_blocks,         // additional blocks with my local info (optional)
                     const diy::Master&         master,
                     diy::mpi::communicator&    comm)
{
    std::vector<int> msg;

    // for all procs in the neighbor info
    for (auto proc_it = shared_neigh_blocks.begin(); proc_it != shared_neigh_blocks.end(); proc_it++)
    {
        int src_proc = proc_it->first;

        // for all messages arriving from this proc (vertex, block, block ...)
        for (auto neigh_it = proc_it->second.begin(); neigh_it != proc_it->second.end(); neigh_it++)
        {
            msg.clear();
            comm.recv(src_proc, 0, msg);

            // sanity check
            if (msg.size() < 2)
            {
                fmt::print(stderr, "Error: recv_neigh_info() message size is {} ints (< 2)\n", msg.size());
                abort();
            }

            // debug
//             fmt::print(stderr, "received from proc {}: [{}]\n", src_proc, fmt::join(msg, ","));

            // match received info with my own info by vertex id (globally unique)
            int vid = msg[0];                       // vertex gid

            // check shared_neigh_blocks for a matching vid
            auto vert_it = proc_it->second.find(vid);
            if (vert_it != proc_it->second.end())
            {
                // debug
//                 fmt::print(stderr, "recv_neigh_info() matched message [{}] with shared_neigh_blocks [{}]\n",
//                         fmt::join(msg, ","), fmt::join(vert_it->second, ","));

                // update the link for my local blocks matching the message
                diy::BlockID block_id;
                for (auto block_it = vert_it->second.begin(); block_it != vert_it->second.end(); block_it++)
                {
                    diy::Link* link = master.link(master.lid(*block_it));

                    for (auto i = 1; i < msg.size(); i++)       // msg[0] is the vertex id, received block gids start at msg[1]
                    {
                        if (link->find(msg[i]) == -1)           // block not linked yet
                        {
                            block_id.gid = msg[i];
                            block_id.proc = src_proc;
                            link->add_neighbor(block_id);
                        }
                    }
                }
            }

            // check extra_neigh_blocks for a matching vid
            if (extra_neigh_blocks.size())
            {
                auto extra_proc_it = extra_neigh_blocks.find(src_proc);
                if (extra_proc_it != extra_neigh_blocks.end())
                {
                    auto extra_vert_it = extra_proc_it->second.find(vid);
                    if (extra_vert_it != extra_proc_it->second.end())
                    {
                        // debug
//                         fmt::print(stderr, "recv_neigh_info() matched message [{}] with extra_neigh_blocks [{}]\n",
//                                 fmt::join(msg, ","), fmt::join(extra_vert_it->second, ","));

                        // update the link for my local blocks matching the message
                        diy::BlockID block_id;
                        for (auto block_it = extra_vert_it->second.begin(); block_it != extra_vert_it->second.end(); block_it++)
                        {
                            diy::Link* link = master.link(master.lid(*block_it));

                            for (auto i = 1; i < msg.size(); i++)       // msg[0] is the vertex id, received block gids start at msg[1]
                            {
                                if (link->find(msg[i]) == -1)           // block not linked yet
                                {
                                    block_id.gid = msg[i];
                                    block_id.proc = src_proc;
                                    link->add_neighbor(block_id);
                                }
                            }
                        }   // for (auto block_it ...
                    }   // if (extra_vert ...
                }   // if (extra_proc_it ...
            }   // if (extra_neigh_blocks.size()
        }   // for all messages arriving from this proc
    }   // for all procs
}

int main(int argc, char**argv)
{
    // initialize MPI
    diy::mpi::environment  env(argc, argv);     // equivalent of MPI_Init(argc, argv)/MPI_Finalize()
    diy::mpi::communicator world;               // equivalent of MPI_COMM_WORLD

    // command line options
    string  input = "tet";                      // input data: "tet" synthetic dense tet mesh, "hex" synthetic dense hex mesh, "file" input file
    int     mesh_size = 10;                     // source mesh size per side
    string  infile;                             // file name in case of file input
    bool    help;

    using namespace opts;
    Options ops;
    ops
        >> Option('h', "help",      help,           "show help")
        >> Option('i', "input",     input,          "input data type 'tet', 'hex', or 'file'")
        >> Option('f', "filename",  infile,         "input file name in case of file input data")
        >> Option('s', "size",      mesh_size,      "mesh size per side in case of tet or hex input data")
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
    std::string read_opts   = "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;DEBUG_IO=0;";
    std::string outfile     = "debug.h5m";
    std::string write_opts  = "PARALLEL=WRITE_PART;DEBUG_IO=0";

    // create moab mesh
    Interface*                      mbi = new Core();                       // moab interface
    ParallelComm*                   pc  = new ParallelComm(mbi, world);     // moab communicator
    EntityHandle                    root;
    ErrorCode                       rval;
    rval = mbi->create_meshset(MESHSET_SET, root); ERR;

    if (input == "tet")
        PrepMesh(1, mesh_size, 0, mbi, pc, root, 1.0, false);
    else if (input == "hex")
        PrepMesh(0, mesh_size, 0, mbi, pc, root, 1.0, false);
    else if (input == "file")
    {
        rval = mbi->load_file(infile.c_str(), &root, read_opts.c_str() ); ERR;
    }
    else
    {
        fmt::print(stderr, "Unrecognized input type\n");
        abort();
    }

    // debug
    // write vtk file, one per process
//     std::string vtk_filename = "debug" + std::to_string(world.rank()) + ".vtk";
//     rval = mbi->write_file(vtk_filename.c_str(), "VTK"); ERR;

    // debug
    // write h5m file, one shared file
//     rval = mbi->write_file("debug.h5m", 0, write_opts.c_str(), &root, 1); ERR;

    // query number of local parts in the parallel moab partition
    Range   parts;
    Tag     part_tag;
    int     nblocks;            // local number of blocks = local number of parts
    int     tot_nblocks;        // total global number of blocks
    rval = mbi->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, part_tag); ERR;
    rval = mbi->get_entities_by_type_and_tag(0, MBENTITYSET, &part_tag, NULL, 1, parts); ERR;
    nblocks = parts.size();

    // sum all local parts to get global total number of blocks
    diy::mpi::all_reduce(world, nblocks, tot_nblocks, std::plus<int>());
    fmt::print(stderr, "number of local blocks nblocks = {} number of global blocks tot_nblocks = {}\n", nblocks, tot_nblocks);

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

    // build map of entities -> parts (blocks)
    EntityPart      entity_part_map;                        // entity global id -> part (block) gid
    build_entity_part_map(parts, part_tag, mbi, pc, entity_part_map);

    // loop over local blocks, filling in their entities, adding blocks to the master, assigner, creating their links
    collect_local_blocks(parts, part_tag, mbi, entity_part_map, master, assigner);

    // get entities owned by my process and shared with other processors
    Range shared_ents;
    rval = pc->get_shared_entities(-1, shared_ents, MBVERTEX); ERR;

    // shared entities owned by my process
    Range my_shared_ents;
    rval = pc->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &my_shared_ents ); ERR;

    // shared entities owned by other processes
    Range other_shared_ents;
    rval = pc->filter_pstatus(shared_ents, PSTATUS_NOT_OWNED, PSTATUS_AND, -1, &other_shared_ents); ERR;

    // collect neighbor info for shared entities owned by my process
    ProcNeighBlocks my_neigh_info;
    collect_neigh_info(my_shared_ents, mbi, pc, entity_part_map, my_neigh_info);

    // collect neighbor info for shared entities owned by other procs
    ProcNeighBlocks other_neigh_info;
    collect_neigh_info(other_shared_ents, mbi, pc, entity_part_map, other_neigh_info);

    // exchange neighbor info with other procs
    send_neigh_info(my_neigh_info, world);
    send_neigh_info(other_neigh_info, world);
    recv_neigh_info(my_neigh_info, other_neigh_info, master, world);
    recv_neigh_info(other_neigh_info, my_neigh_info, master, world);

    // debug: print the link after remote exchange
    for (auto i = 0; i < master.size(); i++)
    {
        Block*      b    = static_cast<Block*>(master.block(i));
        diy::Link*  link = master.link(i);
        fmt::print(stderr, "Final link for block gid {} has size {}:\n", b->gid, link->size());
        for (auto i = 0; i < link->size(); i++)
            fmt::print(stderr, "[gid, proc] = [{}, {}]\n", link->target(i).gid, link->target(i).proc);
    }

    return 0;
}

