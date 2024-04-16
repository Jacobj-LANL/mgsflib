function [mesh,tetmesh] = mgsfmesh_extract_tets(mesh,degree,distribution,nset,fset)

%#codegen -args {SfeMesh, int32(0), int32(0), int32(0), int32(0)}

if coder.target('MATLAB') && exist(['meshing_gateway.' mexext], 'file')
    [mesh,tetmesh] = meshing_gateway('mgsfmesh_extract_tets', mesh,degree,distribution,nset,fset);
    return
end

tetnodes2edges = int32([0, 1, 3, 4; 1, 0, 2, 5; 3, 2, 0, 6; 4, 5, 6, 0]);

etype_hi = encode_elemtype(SFE_SHAPE_TET, degree, distribution);

nv = sfemesh_nnodes(mesh);
bndmask = false(nv, 1);
for i = 1:sfemesh_nnodesets(mesh, nset); bndmask(mesh.nodesets(nset).nids(i)) = true; end

nfacets = sfemesh_nfacetsets(mesh, fset);

% determining which edges are on boundary
mesh = sfemesh_append_edgeset(mesh, 'bndedges');
mesh = sfemesh_resize_edgesetdata(mesh, int32(1), 3 * nfacets);

lnodes = zeros(1, 4);
netypes = sfemesh_nelemtypes(mesh);
nedges = int32(0);
etable_tet = int32(1);
for e = 1:netypes
    etype = mesh.elemtables(e).etype;
    shape = obtain_elemshape(etype);
    if shape == SFE_SHAPE_TET
        etable_tet = e;
        dist = obtain_elemnodepos(etype);
        for ii = 1:sfemesh_nelems(mesh, e)
            nfound = int32(0);
            for jj = int32(1):4
                v = mesh.elemtables(e).conn(ii, jj);
                if bndmask(v)
                    nfound = nfound + 1;
                    lnodes(nfound) = jj;
                end
            end

            if nfound == 2
                edge = tetnodes2edges(lnodes(1), lnodes(2));
                nedges = nedges + 1;
                [~, geid] = sfemesh_leid2teid(e, ii, mesh);
                mesh.edgesets(1).eids(nedges) = geid;
                mesh.edgesets(1).edges(nedges) = edge;
            end
        end
    end
end

mesh = sfemesh_resize_edgesetdata(mesh, int32(1), nedges);
sfes = SfeObject(SFE_TET_4);
natcoords_ = obtain_natcoords(etype_hi);
nnodes = obtain_elemnnodes(etype_hi);
ntets = sfemesh_nelems(mesh,etable_tet);

tetmesh = sfemesh_create(int32(3),nnodes*(nedges + nfacets));
tetmesh = sfemesh_append_econntable(tetmesh, etype_hi);
tetmesh = sfemesh_resize_econn(tetmesh,int32(1),nedges+nfacets);
tetmesh = sfemesh_append_elemset(tetmesh);
tetmesh = sfemesh_resize_elemseteids(tetmesh, int32(1), ntets);

ne = int32(0);
nv = int32(0);
for ii = 1:nfacets
    geid = mesh.facetsets(fset).eids(ii);
    [etable, eid] = sfemesh_teid2leid(mesh, geid);
    ne = ne+1;
    tetmesh.elemsets(1).eids(eid) = ne;
    for jj = 1:nnodes
        nv = nv+1;
        tetmesh.elemtables(1).conn(ne,jj) = nv;
    end
end

for ii = 1:nedges
    geid = mesh.edgesets(1).eids(ii);
    [etable, eid] = sfemesh_teid2leid(mesh, geid);
    ne = ne+1;
    tetmesh.elemsets(1).eids(eid) = ne;
    for jj = 1:nnodes
        nv = nv+1;
        tetmesh.elemtables(1).conn(ne,jj) = nv;
    end
end


end