function mesh = mgsfmesh_improve_tets_stellar_mixed(mesh)
nv = int32(size(mesh.coords,1));
tetmesh = sfemesh_create(int32(3));
tetmesh = sfemesh_append_econntable(tetmesh, SFE_TET_4);
mesh_buff = sfemesh_create(int32(3));

netypes = sfemesh_nelemtypes(mesh);
for e = int32(1):netypes
    % get necessary element data
    etype = mesh.elemtables(e).etype;
    shape = obtain_elemshape(etype);
    degree = obtain_elemdegree(etype);
    nelems = sfemesh_nelems(mesh, e);

    if shape == SFE_SHAPE_TET && degree == 1
        % isolating tetmesh and running through stellar
        tetmesh = sfemesh_resize_econn(tetmesh, int32(1), nelems);
        tetmesh.elemtables(1).conn = mesh.elemtables(e).conn;
        tetmesh = sfemesh_resize_coords(tetmesh,nv);
        tetmesh.coords = mesh.coords;
        tetmesh = remove_unused_nodes(tetmesh);
        bndtags = sfemesh_determine_bndnodes(tetmesh);
        [tetmesh.elemtables(1).conn, tetmesh.coords] = call_Stellar_fileio(tetmesh.elemtables(1).conn, tetmesh.coords, bndtags);
        tetmesh = remove_unused_nodes(tetmesh);
        bndtags = sfemesh_determine_bndnodes(tetmesh);
        bndnids = int32(find(bndtags));
        
        % isolating the rest of the mesh and finding boundary of pyramids
        mesh_buff.coords = mesh.coords;
        mesh_buff.elemtables = mesh.elemtables(setdiff(1:netypes, e));
        mesh_buff = remove_unused_nodes(mesh_buff);
        bndtags2 = false(size(mesh_buff.coords,1),1);
        for ee = int32(1):sfemesh_nelemtypes(mesh_buff)
            etype2 = mesh_buff.elemtables(ee).etype;
            shape2 = obtain_elemshape(etype2);
            if shape2 == SFE_SHAPE_PYRA
                istart = mesh_buff.elemtables(ee).istart;
                for ii = 1:sfemesh_nelems(mesh_buff,ee)
                    eid = istart+ii-1;
                    for jj = 2:5
                        if ~mesh_buff.sibhfs(eid,jj)
                            [~, lids] = obtain_facets(etype2, int8(jj));
                            nn = cast(size(lids, 1), 'int32');
                            for k = 1:nn
                                kk = cast(lids(k), 'int32');
                                bndtags2(mesh_buff.elemtables(ee).conn(ii,kk)) = true;
                            end
                        end
                    end
                end
            end
        end
        bndnids2 = int32(find(bndtags2));
        
        % matching tet and pyramid boundaries
        nv2 = int32(size(mesh_buff.coords,1));
        pyr2tet = zeros(nv2,1,'int32');
        tet2pyr = zeros(size(tetmesh.coords,1),1,'int32');
        for ii = 1:length(bndnids2)
            v2 = bndnids2(ii);
            for jj = 1:length(bndnids)
                v1 = bndnids(jj);
                if norm(tetmesh.coords(v1,1:3)-mesh_buff.coords(v2,1:3)) < 1e-8
                    pyr2tet(v2) = v1;
                    tet2pyr(v1) = v2;
                    break;
                end
            end
        end

        for ii = 1:sfemesh_nelems(tetmesh,1)
            for jj = 1:4
                v = tetmesh.elemtables(1).conn(ii,jj);
                if tet2pyr(v)
                    tetmesh.elemtables(1).conn(ii,jj) = tet2pyr(v);
                else
                    tetmesh.elemtables(1).conn(ii,jj) = v + nv2;
                end
            end
        end

        mesh.elemtables(setdiff(1:netypes, e)) = mesh_buff.elemtables;
        mesh.elemtables(e) = tetmesh.elemtables(1);
        mesh.coords = [mesh_buff.coords; tetmesh.coords];
        mesh = remove_unused_nodes(mesh);

        

    end
end

end