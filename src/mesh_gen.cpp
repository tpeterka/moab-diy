#include "prod-con.hpp"

using namespace std;
using namespace moab;

// generate a regular structured hex mesh
void hex_mesh_gen(int *mesh_size,                       // mesh size (i,j,k) number of vertices in each dim
        Interface *mbint,                               // moab interface instance
        EntityHandle *mesh_set,                         // moab mesh set
        ParallelComm *mbpc,                             // moab parallel communicator
        diy::RegularDecomposer<Bounds>& decomp,         // diy decomposition
        diy::RoundRobinAssigner& assign)                // diy assignment
{
    create_hexes_and_verts(mesh_size, mbint, mesh_set, decomp, assign, mbpc);
    resolve_and_exchange(mbint, mesh_set, mbpc);
}

// generate a regular structured tet mesh
void tet_mesh_gen(int *mesh_size,                        // mesh size (i,j,k) number of vertices in each dim
        Interface *mbint,                                // moab interface instance
        EntityHandle *mesh_set,                          // moab mesh set
        ParallelComm *mbpc,                              // moab parallel communicator
        diy::RegularDecomposer<Bounds>& decomp,          // diy decomposition
        diy::RoundRobinAssigner& assign)                 // diy assignment
{
    create_tets_and_verts(mesh_size, mbint, mesh_set, decomp, assign, mbpc);
    resolve_and_exchange(mbint, mesh_set, mbpc);
}

// create hex cells and vertices
void create_hexes_and_verts(int *mesh_size,     // mesh size (i,j,k) number of vertices in each dim
        Interface *mbint,                       // moab interface instance
        EntityHandle *mesh_set,                 // moab mesh set
        diy::RegularDecomposer<Bounds>& decomp, // diy decomposition
        diy::RoundRobinAssigner& assign,        // diy assignment
        ParallelComm *mbpc)                     // moab communicator
{
    Core *mbcore = dynamic_cast<Core*>(mbint);
    EntityHandle vs, cs;
    ErrorCode rval;
    EntityHandle handle;
    vector<int> local_gids;
    int rank;
    MPI_Comm_rank(mbpc->comm(), &rank);
    assign.local_gids(rank, local_gids);
    int nblocks = local_gids.size();
    std::vector<Bounds> bounds;
    for (auto i = 0; i < nblocks; i++)
        bounds.push_back(3);
    int b; // local block number

    // get the read interface from moab
    ReadUtilIface *iface;
    rval = mbint->query_interface(iface); ERR;

    // block bounds
    for (int b = 0; b < nblocks; b++)
        decomp.fill_bounds(bounds[b], local_gids[b], false);

    // todo: multiple blocks per process not really supported yet in this example
    b = 0;

    // the following method is based on the example in
    // moab/examples/old/FileRead.cpp, using the ReadUtilIface class

    // allocate a block of vertex handles and store xyz’s into them
    // returns a starting handle for the node sequence
    vector<double*> arrays;
    EntityHandle startv;
    int num_verts =
        (bounds[b].max[0] - bounds[b].min[0] + 1) *
        (bounds[b].max[1] - bounds[b].min[1] + 1) *
        (bounds[b].max[2] - bounds[b].min[2] + 1);
    rval = iface->get_node_coords(3, num_verts, 0, startv, arrays); ERR;

    // populate vertex arrays
    // vertices normalized to be in the range [0.0 - 1.0]
    int n = 0;
    for (int k = bounds[b].min[2]; k <= bounds[b].max[2]; k++)
    {
        for (int j = bounds[b].min[1]; j <= bounds[b].max[1]; j++)
        {
            for (int i = bounds[b].min[0]; i <= bounds[b].max[0]; i++)
            {
                arrays[0][n] = double(i) / (mesh_size[0] - 1);
                arrays[1][n] = double(j) / (mesh_size[1] - 1);
                arrays[2][n] = double(k) / (mesh_size[2] - 1);

                // debug
                // 	fprintf(stderr, "[i,j,k] = [%d %d %d] vert[%d] = [%.2lf %.2lf %.2lf]\n",
                // 		i, j, k, n, arrays[0][n], arrays[1][n], arrays[2][n]);

                n++;
            }
        }
    }

    // allocate connectivity array
    EntityHandle startc;           // handle for start of cells
    EntityHandle *starth;          // handle for start of connectivity
    int num_hexes =
        (bounds[b].max[0] - bounds[b].min[0]) *
        (bounds[b].max[1] - bounds[b].min[1]) *
        (bounds[b].max[2] - bounds[b].min[2]);
    rval = iface->get_element_connect(num_hexes, 8, MBHEX, 0, startc, starth);

    // populate the connectivity array
    n = 0;
    int m = 0;
    for (int k = bounds[b].min[2]; k <= bounds[b].max[2]; k++)
    {
        for (int j = bounds[b].min[1]; j <= bounds[b].max[1]; j++)
        {
            for (int i = bounds[b].min[0]; i <= bounds[b].max[0]; i++)
            {
                if (i < bounds[b].max[0] && j < bounds[b].max[1] && k < bounds[b].max[2])
                {
                    int A, B, C, D, E, F, G, H; // hex verts according to my diagram
                    D = n;
                    C = D + 1;
                    H = D + bounds[b].max[0] - bounds[b].min[0] + 1;
                    G = H + 1;
                    A = D + (bounds[b].max[0] - bounds[b].min[0] + 1) *
                        (bounds[b].max[1] - bounds[b].min[1] + 1);
                    B = A + 1;
                    E = A + bounds[b].max[0] - bounds[b].min[0] + 1;
                    F = E + 1;

                    // hex ABCDEFGH
                    starth[m++] = startv + A;
                    starth[m++] = startv + B;
                    starth[m++] = startv + C;
                    starth[m++] = startv + D;
                    starth[m++] = startv + E;
                    starth[m++] = startv + F;
                    starth[m++] = startv + G;
                    starth[m++] = startv + H;
                }
                n++;
            }
        }
    }

    // add vertices and cells to the mesh set
    Range vRange(startv, startv + num_verts - 1);      // vertex range
    Range cRange(startc, startc + num_hexes - 1);      // cell range
    rval = mbint->add_entities(*mesh_set, vRange); ERR;
    rval = mbint->add_entities(*mesh_set, cRange); ERR;

    // check that long is indeed 8 bytes on this machine
    assert(sizeof(long) == 8);

    // set global ids
    long gid;
    Tag global_id_tag;
    rval = mbint->tag_get_handle("HANDLEID", sizeof(long), MB_TYPE_OPAQUE,
            global_id_tag, MB_TAG_CREAT|MB_TAG_DENSE); ERR;

    // gids for vertices, starting at 1 by moab convention
    handle = startv;
    for (int k = bounds[b].min[2]; k < bounds[b].max[2] + 1; k++)
    {
        for (int j = bounds[b].min[1]; j < bounds[b].max[1] + 1; j++)
        {
            for (int i = bounds[b].min[0]; i < bounds[b].max[0] + 1; i++)
            {
                gid = (long)1 + (long)i + (long)j * (mesh_size[0]) +
                    (long)k * (mesh_size[0]) * (mesh_size[1]);
//                 fprintf(stderr, "i,j,k = [%d %d %d] gid = %ld\n", i, j, k, gid);
                rval = mbint->tag_set_data(global_id_tag, &handle, 1, &gid); ERR;
                handle++;
            }
        }
    }

    // gids for cells, starting at 1 by moab convention
    handle = startc;
    for (int k = bounds[b].min[2]; k < bounds[b].max[2]; k++)
    {
        for (int j = bounds[b].min[1]; j < bounds[b].max[1]; j++)
        {
            for (int i = bounds[b].min[0]; i < bounds[b].max[0]; i++)
            {
                gid = (long)1 + (long)i + (long)j * (mesh_size[0] - 1) +
                    (long)k * (mesh_size[0] - 1) * (mesh_size[1] - 1);
//                 fprintf(stderr, "i,j,k = [%d %d %d] gid = %ld\n", i, j, k, gid);
                rval = mbint->tag_set_data(global_id_tag, &handle, 1, &gid); ERR;
                handle++;
            }
        }
    }

    // Create a part set to add entities in the current block
    // Technically we could have an array of these to represent multiple
    // blocks in the current process - but we will simplify this assumption
    EntityHandle partset;
    rval = mbint->create_meshset(MESHSET_SET, partset); ERR;
    rval = mbint->add_entities(partset, cRange); ERR;

    Tag parttag;
    int dumid=-1, proc_id=mbpc->rank();
    rval = mbint->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, parttag, MB_TAG_CREAT | MB_TAG_SPARSE, &dumid); ERR;
    rval = mbint->tag_set_data(parttag, &partset, 1, &proc_id); ERR;
    rval = mbint->add_entities(*mesh_set, &partset, 1); ERR;

    Range psets;
    psets.insert(partset);
    mbpc->partition_sets() = psets;

    // update adjacencies (needed by moab)
    rval = iface->update_adjacencies(startc, num_hexes, 8, starth); ERR;

    // cleanup
    rval = mbint->release_interface(iface); ERR;
}

// create tet cells and vertices
void create_tets_and_verts(int *mesh_size,      // mesh size (i,j,k) number of vertices in each dim
        Interface *mbint,                       // moab interface instance
        EntityHandle *mesh_set,                 // moab parallel communicator
        diy::RegularDecomposer<Bounds>& decomp, // diy decomposition
        diy::RoundRobinAssigner& assign,        // diy assignment
        ParallelComm *mbpc)                     // moab communicator
{
    Core *mbcore = dynamic_cast<Core*>(mbint);
    EntityHandle vs, cs;
    ErrorCode rval;
    EntityHandle handle;
    vector<int> local_gids;
    int rank;
    MPI_Comm_rank(mbpc->comm(), &rank);
    assign.local_gids(rank, local_gids);
    int nblocks = local_gids.size();
    std::vector<Bounds> bounds;
    for (auto i = 0; i < nblocks; i++)
        bounds.push_back(3);
    int b;                                   // local block number

    // get the read interface from moab
    ReadUtilIface *iface;
    rval = mbint->query_interface(iface); ERR;

    // block bounds
    for (int b = 0; b < nblocks; b++)
        decomp.fill_bounds(bounds[b], local_gids[b], false);

    // debug
    //     fprintf(stderr, "nblocks = %d\n", nblocks);
    //     fprintf(stderr, "bounds min = [%d %d %d] max = [%d %d %d]\n",
    //             bounds[0].min[0], bounds[0].min[1], bounds[0].min[2],
    //             bounds[0].max[0], bounds[0].max[1], bounds[0].max[2]);

    // todo: multiple blocks per process not really supported yet in this example
    b = 0;

    // the following method is based on the example in
    // moab/examples/old/FileRead.cpp, using the ReadUtilIface class

    // allocate a block of vertex handles and store xyz’s into them
    // returns a starting handle for the node sequence
    vector<double*> arrays;
    EntityHandle startv;
    int num_verts =
        (bounds[b].max[0] - bounds[b].min[0] + 1) *
        (bounds[b].max[1] - bounds[b].min[1] + 1) *
        (bounds[b].max[2] - bounds[b].min[2] + 1);
    rval = iface->get_node_coords(3, num_verts, 0, startv, arrays); ERR;

    // populate vertex arrays
    // vertices normalized to be in the range [0.0 - 1.0]
    int n = 0;
    for(int k = bounds[b].min[2]; k <= bounds[b].max[2]; k++)
    {
        for(int j = bounds[b].min[1]; j <= bounds[b].max[1]; j++)
        {
            for(int i = bounds[b].min[0]; i <= bounds[b].max[0]; i++)
            {
                arrays[0][n] = double(i) / (mesh_size[0] - 1);
                arrays[1][n] = double(j) / (mesh_size[1] - 1);
                arrays[2][n] = double(k) / (mesh_size[2] - 1);

                // debug
                //                 fprintf(stderr, "vert[%d] = [%.2lf %.2lf %.2lf]\n",
                //                         n, arrays[0][n], arrays[1][n], arrays[2][n]);
                n++;
            }
        }
    }

    // allocate connectivity array
    EntityHandle startc;                          // handle for start of cells
    EntityHandle *starth;                         // handle for start of connectivity
    int num_tets = 6 *                            // each hex cell will be converted to 6 tets
        (bounds[b].max[0] - bounds[b].min[0]) *
        (bounds[b].max[1] - bounds[b].min[1]) *
        (bounds[b].max[2] - bounds[b].min[2]);
    rval = iface->get_element_connect(num_tets, 4, MBTET, 0, startc, starth); ERR;

    // populate the connectivity array
    n = 0;
    int m = 0;
    for (int k = bounds[b].min[2]; k <= bounds[b].max[2]; k++)
    {
        for (int j = bounds[b].min[1]; j <= bounds[b].max[1]; j++)
        {
            for (int i = bounds[b].min[0]; i <= bounds[b].max[0]; i++)
            {
                if (i < bounds[b].max[0] && j < bounds[b].max[1] &&
                        k < bounds[b].max[2])
                {
                    int A, B, C, D, E, F, G, H;   // hex verts according to my diagram
                    D = n;
                    C = D + 1;
                    H = D + bounds[b].max[0] - bounds[b].min[0] + 1;
                    G = H + 1;
                    A = D + (bounds[b].max[0] - bounds[b].min[0] + 1) *
                        (bounds[b].max[1] - bounds[b].min[1] + 1);
                    B = A + 1;
                    E = A + bounds[b].max[0] - bounds[b].min[0] + 1;
                    F = E + 1;

                    // tet EDHG
                    starth[m++] = startv + E;
                    starth[m++] = startv + D;
                    starth[m++] = startv + H;
                    starth[m++] = startv + G;

                    // tet ABCF
                    starth[m++] = startv + A;
                    starth[m++] = startv + B;
                    starth[m++] = startv + C;
                    starth[m++] = startv + F;

                    // tet ADEF
                    starth[m++] = startv + A;
                    starth[m++] = startv + D;
                    starth[m++] = startv + E;
                    starth[m++] = startv + F;

                    // tet CGDF
                    starth[m++] = startv + C;
                    starth[m++] = startv + G;
                    starth[m++] = startv + D;
                    starth[m++] = startv + F;

                    // tet ACDF
                    starth[m++] = startv + A;
                    starth[m++] = startv + C;
                    starth[m++] = startv + D;
                    starth[m++] = startv + F;

                    // tet DGEF
                    starth[m++] = startv + D;
                    starth[m++] = startv + G;
                    starth[m++] = startv + E;
                    starth[m++] = startv + F;
                }
                n++;
            }
        }
    }

    // add vertices and cells to the mesh set
    Range vRange(startv, startv + num_verts - 1);      // vertex range
    Range cRange(startc, startc + num_tets - 1);       // cell range
    rval = mbint->add_entities(*mesh_set, vRange); ERR;
    rval = mbint->add_entities(*mesh_set, cRange); ERR;

    // set global ids
    long gid;
    Tag global_id_tag;
    rval = mbint->tag_get_handle("HANDLEID", sizeof(long), MB_TYPE_OPAQUE,
            global_id_tag, MB_TAG_CREAT|MB_TAG_DENSE); ERR;

    // gids for vertices, starting at 1 by moab convention
    handle = startv;
    for (int k = bounds[b].min[2]; k < bounds[b].max[2] + 1; k++)
    {
        for (int j = bounds[b].min[1]; j < bounds[b].max[1] + 1; j++)
        {
            for (int i = bounds[b].min[0]; i < bounds[b].max[0] + 1; i++)
            {
                gid = (long)1 + (long)i + (long)j * (mesh_size[0]) +
                    (long)k * (mesh_size[0]) * (mesh_size[1]);
                // debug
                //                 fprintf(stderr, "i,j,k = [%d %d %d] gid = %ld\n", i, j, k, gid);
                rval = mbint->tag_set_data(global_id_tag, &handle, 1, &gid); ERR;
                handle++;
            }
        }
    }

    // gids for cells, starting at 1 by moab convention
    handle = startc;
    for (int k = bounds[b].min[2]; k < bounds[b].max[2]; k++)
    {
        for (int j = bounds[b].min[1]; j < bounds[b].max[1]; j++)
        {
            for (int i = bounds[b].min[0]; i < bounds[b].max[0]; i++)
            {
                for (int t = 0; t < 6; t++)            // 6 tets per grid space
                {
                    gid = (long)1 + (long)t +  (long)i * 6 + (long)j * 6 * (mesh_size[0] - 1) +
                        (long)k * 6 * (mesh_size[0] - 1) * (mesh_size[1] - 1);
                    // 	 fprintf(stderr, "t,i,j,k = [%d %d %d %d] gid = %ld\n", t, i, j, k, gid);
                    rval = mbint->tag_set_data(global_id_tag, &handle, 1, &gid); ERR;
                    handle++;
                }
            }
        }
    }

    // update adjacencies (needed by moab)
    rval = iface->update_adjacencies(startc, num_tets, 4, starth); ERR;

    // Create a part set to add entities in the current block
    // Technically we could have an array of these to represent multiple
    // blocks in the current process - but we will simplify this assumption
    EntityHandle partset;
    rval = mbint->create_meshset(MESHSET_SET, partset); ERR;
    rval = mbint->add_entities(partset, cRange); ERR;

    Tag parttag;
    int dumid=-1, proc_id=mbpc->rank();
    rval = mbint->tag_get_handle("PARALLEL_PARTITION", 1, MB_TYPE_INTEGER, parttag, MB_TAG_CREAT | MB_TAG_SPARSE, &dumid); ERR;
    rval = mbint->tag_set_data(parttag, &partset, 1, &proc_id); ERR;
    rval = mbint->add_entities(*mesh_set, &partset, 1); ERR;

    Range psets;
    psets.insert(partset);
    mbpc->partition_sets() = psets;

    // cleanup
    rval = mbint->release_interface(iface); ERR;

}

// resolve shared entities
void resolve_and_exchange(Interface *mbint,       // mbint: moab interface instance
        EntityHandle *mesh_set, // mesh_set: moab mesh set
        ParallelComm *mbpc)     // mbpc: moab parallel communicator
{
    ErrorCode rval;

    mbpc->partition_sets().insert(*mesh_set);
    Tag global_id_tag;
    rval = mbint->tag_get_handle("HANDLEID", sizeof(long), MB_TYPE_OPAQUE, global_id_tag, MB_TAG_DENSE); ERR;
    rval = mbpc->resolve_shared_ents(*mesh_set, -1, -1, &global_id_tag); ERR;
}

// return a value for the field position (simple magnitude)
double PhysField(double x,
        double y,
        double z,
        double factor)
{
    return factor * sqrt(x * x + y * y + z * z);
}

// add a value to each element in the field
void PutElementField(Interface *mbi,
        EntityHandle eh,
        const char *tagname,
        double factor)
{
    Range elems;
    ErrorCode rval;
    const double defVal = 0.;
    Tag fieldTag;

    rval = mbi->get_entities_by_dimension(eh, 3, elems); ERR;
    rval = mbi->tag_get_handle(tagname, 1, MB_TYPE_DOUBLE, fieldTag,
            MB_TAG_DENSE|MB_TAG_CREAT, &defVal); ERR;

    for(int i = 0; i< (int)elems.size(); i++)
    {
        EntityHandle elem = elems[i];
        double pos[3];
        rval = mbi->get_coords(&elem, 1, pos); ERR;
        double field_value =  PhysField(pos[0], pos[1], pos[2], factor);
        rval = mbi->tag_set_data(fieldTag, &elem, 1, &field_value); ERR;

        // debug
        //     fprintf(stderr, "setting element %d: f(%.3lf, %.3lf, %.3lf) = %.3lf\n",
        // 	    i, x, y, z, fieldValue);
    }
}

// gets the element field
void GetElementField(Interface *mbi,
        EntityHandle eh,
        const char *tagname,
        double factor,
        MPI_Comm comm,
        bool debug)
{
    Range elems;
    double *field_values;

    mbi->get_entities_by_dimension(eh, 3, elems);

    Tag field_tag;
    int count;
    mbi->tag_get_handle(tagname, 1, MB_TYPE_DOUBLE, field_tag);
    mbi->tag_iterate(field_tag, elems.begin(), elems.end(), count,
            (void *&)field_values);
    assert(count == (int)elems.size());

    for(int i = 0; i < (int)elems.size(); i++)
    {
        EntityHandle elem = elems[i];
        double pos[3];
        mbi->get_coords(&elem, 1, pos);
        double ref_value = PhysField(pos[0], pos[1], pos[2], factor);

        // debug
        if (debug)
            fprintf(stderr, "element %d: f(%.3lf, %.3lf, %.3lf) = %.3lf "
                    "ref_value = %.3lf\n",
                    i, pos[0], pos[1], pos[2], field_values[i], ref_value);
    }
}

// add a value to each vertex in the field
void PutVertexField(Interface *mbi,
        EntityHandle eh,
        const char *tagname,
        double factor)
{
    Range verts;

    mbi->get_entities_by_type(eh, MBVERTEX, verts);

    const double defVal = 0.;
    Tag fieldTag;
    mbi->tag_get_handle(tagname, 1, MB_TYPE_DOUBLE, fieldTag,
            MB_TAG_DENSE|MB_TAG_CREAT, &defVal);

    int numVerts = verts.size();

    for(int i = 0; i < numVerts; i++)
    {
        EntityHandle vert = verts[i];
        double pos[3];
        mbi->get_coords(&vert, 1, pos);
        double fieldValue =  PhysField(pos[0], pos[1], pos[2], factor);
        mbi->tag_set_data(fieldTag, &vert, 1, &fieldValue);

        // debug
        //     fprintf(stderr, "setting element %d: f(%.3lf, %.3lf, %.3lf) = %.3lf\n",
        // 	    i, pos[0], pos[1], pos[2], fieldValue);
    }
}

// gets the vertex field
void GetVertexField(Interface *mbi,
        EntityHandle eh,
        const char *tagname,
        double factor,
        MPI_Comm comm,
        bool debug)
{
    Range verts;
    double *field_values;

    mbi->get_entities_by_type(eh, MBVERTEX, verts);

    Tag field_tag;
    int count;
    mbi->tag_get_handle(tagname, 1, MB_TYPE_DOUBLE, field_tag);
    mbi->tag_iterate(field_tag, verts.begin(), verts.end(), count, (void *&)field_values);
    assert(count == (int)verts.size());

    for(int i = 0; i < (int)verts.size(); i++)
    {
        EntityHandle vert = verts[i];
        double pos[3];
        mbi->get_coords(&vert, 1, pos);
        double ref_value = PhysField(pos[0], pos[1], pos[2], factor);

        // debug
        if (debug)
            fprintf(stderr, "vertex %d: f(%.3lf, %.3lf, %.3lf) = %.3lf "
                    "ref_value = %.3lf\n",
                    i, pos[0], pos[1], pos[2], field_values[i], ref_value);
    }
}
// prints mesh statistics (for debugging)
void PrintMeshStats(Interface *mbint,        // moab interface
        EntityHandle *mesh_set,  // moab mesh set
        ParallelComm *mbpc)      // moab parallel communicator
{
    Range range;
    ErrorCode rval;
    static int mesh_num = 0;                 // counts how many time this function is called
    float loc_verts = 0.0;                   // num local verts (fractional for shared verts)
    float glo_verts;                         // global number of verts
    float loc_share_verts = 0.0;             // num local shared verts (fractional)
    float glo_share_verts;                   // global number of shared verts
    int loc_cells, glo_cells;                // local and global number of cells (no sharing)
    int rank;

    MPI_Comm_rank(mbpc->comm(), &rank);

    // get local quantities
    range.clear();
    rval = mbint->get_entities_by_dimension(*mesh_set, 0, range); ERR;

    // compute fractional contribution of shared vertices attributed to this proc
    int ps[MAX_SHARING_PROCS];               // sharing procs for a vert
    EntityHandle hs[MAX_SHARING_PROCS];      // handles of shared vert on sharing procs
    int num_ps = 0;                          // number of sharing procs, returned by moab
    unsigned char pstat;                     // pstatus, returned by moab
    for (Range::iterator verts_it = range.begin(); verts_it != range.end();
            verts_it++)
    {
        rval = mbpc->get_sharing_data(*verts_it, ps, hs, pstat, num_ps); ERR;
        if (num_ps == 0)
            loc_verts++;
        else if (num_ps == 1)                // when 2 procs, moab lists only one (the other)
        {
            loc_verts += 0.5;
            loc_share_verts += 0.5;
        }
        else
        {
            loc_verts += (1.0 / (float)num_ps);
            loc_share_verts += (1.0 / (float)num_ps);
        }
    }

    range.clear();
    rval = mbint->get_entities_by_dimension(*mesh_set, 3, range); ERR;
    loc_cells = (int)range.size();

    // add totals for global quantities
    MPI_Reduce(&loc_verts, &glo_verts, 1, MPI_FLOAT, MPI_SUM, 0,
            mbpc->comm());
    MPI_Reduce(&loc_share_verts, &glo_share_verts, 1, MPI_FLOAT, MPI_SUM, 0,
            mbpc->comm());
    MPI_Reduce(&loc_cells, &glo_cells, 1, MPI_INT, MPI_SUM, 0, mbpc->comm());

    // report results
    if (rank == 0)
    {
        fmt::print(stderr, "----------------- Mesh {} statistics -----------------\n", mesh_num);
        fmt::print(stderr, "Total number of verts = {} of which {} "
                "are shared\n", glo_verts, glo_share_verts);
        fmt::print(stderr, "Total number of cells = {}\n", glo_cells);
        fmt::print(stderr, "------------------------------------------------------\n");
    }

    mesh_num = (mesh_num + 1) % 2;
}

// prepares the mesh by decomposing source domain and creating mesh in situ
void PrepMesh(int src_type,
              int src_size,
              int slab,
              Interface* mbi,
              ParallelComm* pc,
              EntityHandle root,
              double factor,
              bool debug)
{
    diy::RegularDecomposer<Bounds>::CoordinateVector    ghost(3, 0);
    diy::RegularDecomposer<Bounds>::DivisionsVector     given(3, 0);
    std::vector<bool>                                   share_face(3, true);
    std::vector<bool>                                   wrap(3, false);

    diy::mpi::communicator comm = pc->comm();

    int src_mesh_size[3]  = {src_size, src_size, src_size};      // source size
    Bounds domain(3);
    domain.min[0] = domain.min[1] = domain.min[2] = 0;
    domain.max[0] = domain.max[1] = domain.max[2] = src_size - 1;

    // decompose domain
    if (slab == 1)
    {
        // blocks are slabs in the z direction
        given[1] = 1;
    }

    // the following calls to asssigners are for 1 block per process
    diy::RoundRobinAssigner         assigner(comm.size(), comm.size());
    diy::RegularDecomposer<Bounds>  decomposer(3,
                                               domain,
                                               assigner.nblocks(),
                                               share_face,
                                               wrap,
                                               ghost,
                                               given);

    // report the number of blocks in each dimension of each mesh
    if (comm.rank() == 0)
        fmt::print(stderr, "Number of blocks in source = [{}]\n", fmt::join(decomposer.divisions, ","));

    // create mesh in situ
    if (debug)
        pc->set_debug_verbosity(5);

    if (src_type == 0)
        hex_mesh_gen(src_mesh_size, mbi, &root, pc, decomposer, assigner);
    else
        tet_mesh_gen(src_mesh_size, mbi, &root, pc, decomposer, assigner);


    // debug: print mesh stats
    PrintMeshStats(mbi, &root, pc);

    // add field to input mesh
    PutVertexField(mbi, root, "vertex_field", factor);
    PutElementField(mbi, root, "element_field", factor);
}


