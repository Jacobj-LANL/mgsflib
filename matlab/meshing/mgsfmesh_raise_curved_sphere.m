function [mesh, tetmesh] = mgsfmesh_raise_curved_sphere(mesh, tetmesh, r)

%#codegen -args {SfeMesh, SfeMesh, 0}
if coder.target('MATLAB') && exist(['meshing_gateway.' mexext], 'file')
    [mesh, tetmesh] = meshing_gateway('mgsfmesh_raise_curved_sphere', mesh, tetmesh, r);
    return
end

geom_degree = obtain_elemdegree(tetmesh.elemtables(1).etype);

nfacets = sfemesh_nfacetsets(mesh, int32(7));
nedges = sfemesh_nedgesets(mesh, int32(1));

% raise facets to be curved
xsloc_ = zeros(84,3);
etype = mesh.elemtables(3).etype;
ndims = int32(size(mesh.coords,2));
for ii = 1:nfacets
    [etable, eid] = sfemesh_teid2leid(mesh, mesh.facetsets(7).eids(ii));
    etype = mesh.elemtables(etable).etype;
    dist = obtain_elemnodepos(etype);
    shape = obtain_elemshape(etype);
    nnodes_hi = obtain_elemnnodes(etype);
    facet = mesh.facetsets(7).facets(ii);
    xsloc_(1:4,1:3) = mesh.coords(mesh.elemtables(etable).conn(eid, 1:4), 1:3);
    sfes = SfeObject(etype);
    etype_hi = etype;
    
    for deg = int32(1):geom_degree - 1
        etype_low = encode_elemtype(shape, deg, int32(deg > 2 && deg < 5 && dist));
        nnodes = obtain_elemnnodes(etype_low);
        etype_hi = encode_elemtype(shape, deg + 1, int32(deg + 1 > 2 && deg + 1 < 5 && dist));
        nnodes_hi = obtain_elemnnodes(etype_hi);
        natcoords = obtain_natcoords3d(etype_hi);

        sfes = sfe_init(sfes, etype_low, xsloc_(1:nnodes,1:3), natcoords);
        xsloc_(1:nnodes_hi,1:3) = sfes.cs_phy;

        [~, lids] = obtain_facets(etype_hi, facet);

        % projecting onto the sphere
        for nn = int32(4):length(lids)
            xsloc_(lids(nn), 1:3) = project_onto_sphere(xsloc_(lids(nn), 1:3), r);
        end
    end

    % inserting data in tetmesh
    ne = tetmesh.elemsets(1).eids(eid);
    for i = 1:nnodes_hi
        for j = 1:ndims
            tetmesh.coords(tetmesh.elemtables(int32(1)).conn(ne, i), j) = xsloc_(i,j);
        end
    end

    % interpolating down to mesh element degree
    sfes = sfe_init(sfes, etype_hi, xsloc_(1:nnodes_hi,1:3), obtain_natcoords3d(etype));
    mesh.coords(mesh.elemtables(etable).conn(eid, 5:end), 1:3) = sfes.cs_phy(5:end, 1:3);
end

% raise edges to be curved
for ii = 1:nedges
    geid = mesh.edgesets(1).eids(ii);
    [etable, eid] = sfemesh_teid2leid(mesh, geid);
    etype = mesh.elemtables(etable).etype;
    dist = obtain_elemnodepos(etype);
    shape = obtain_elemshape(etype);
    nnodes_hi = obtain_elemnnodes(etype);
    edge = mesh.edgesets(1).edges(ii);
    xsloc_(1:4,1:3) = mesh.coords(mesh.elemtables(etable).conn(eid, 1:4), 1:3);
    sfes = SfeObject(etype);
    etype_hi = etype;

    for deg = int32(1):geom_degree - 1
        etype_low = encode_elemtype(shape, deg, int32(deg > 2 && deg < 5 && dist));
        nnodes = obtain_elemnnodes(etype_low);
        etype_hi = encode_elemtype(shape, deg + 1, int32(deg + 1 > 2 && deg + 1 < 5 && dist));
        nnodes_hi = obtain_elemnnodes(etype_hi);
        natcoords = obtain_natcoords3d(etype_hi);

        sfes = sfe_init(sfes, etype_low, xsloc_(1:nnodes,1:3), natcoords);
        xsloc_(1:nnodes_hi,1:3) = sfes.cs_phy;

        [~, lids] = obtain_edges_tet(etype_hi, edge);

        % projecting onto the sphere
        for nn = int32(3):length(lids)
            xsloc_(lids(nn), 1:3) = project_onto_sphere(xsloc_(lids(nn), 1:3), r);
        end
    end

    % inserting data in tetmesh
    ne = tetmesh.elemsets(1).eids(eid);
    for i = 1:nnodes_hi
    for j = 1:ndims; tetmesh.coords(tetmesh.elemtables(1).conn(ne, i), j) = xsloc_(i,j); end
    end

    % interpolating down to mesh element degree
    sfes = sfe_init(sfes, etype_hi, xsloc_(1:nnodes_hi,1:3), obtain_natcoords3d(etype));
    mesh.coords(mesh.elemtables(etable).conn(eid, 5:end), 1:3) = sfes.cs_phy(5:end, 1:3);
end

end

function xs_proj = project_onto_sphere(xs, r)
coder.inline('always');
a = r / sqrt((xs(1) - 0.5)^2 + (xs(2) - 0.5)^2 + (xs(3) - 0.5)^2);
xs_proj = [0, 0, 0];
xs_proj(1) = a * (xs(1) - 0.5) + 0.5;
xs_proj(2) = a * (xs(2) - 0.5) + 0.5;
xs_proj(3) = a * (xs(3) - 0.5) + 0.5;
end
