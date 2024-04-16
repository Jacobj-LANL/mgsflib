function [mesh, aes_mesh, aes2mesh,mesh2aes] = mgsfmesh_extract_mesh(mesh, nlayers)
% mgsfmesh_extract_mesh - extracts 3d aes-fem mesh from standard fem mesh
%
%  [mesh, aes_mesh, aes2mesh] = mgsfmesh_extract_mesh(mesh, nlayers)
%
% Parameters
% ----------
%   mesh:      A WlsMesh or SfeMesh object passed by reference
%   nlayers:    Number of layers of linear quadrilaterals extracted
%
% Returns
% -------
%   mesh:    SfeMesh object passed by reference
%   aes_mesh: SfeMesh object containing only linear tetrahedra
%
% Notes
% -----
% The aes-fem mesh containg decomposed linear tets from input mesh
% plus the decomposed and split pyramids and hexahedrals share a face with
% a pyramid

%#codegen -args {SfeMesh, int32(0)}


if coder.target("MATLAB")
if nargin < 1
    [mesh,tetmesh] = mgsfmesh_shape_mixed('sphere', 20, 3, 'Equidistant', 0.1, false);
    nlayers = int32(3);
end
end


if coder.target('MATLAB') && exist(['meshing_gateway.' mexext], 'file')
    [mesh,aes_mesh, aes2mesh,mesh2aes] = meshing_gateway('mgsfmesh_extract_mesh', mesh, nlayers);
    return
end

m2cAssert(mesh.topo_ndims == 3);
degree = obtain_elemdegree(mesh.elemtables(1).etype); % assuming all elements have same degree

nlayersbig = int32(max(ceil(double(nlayers)/double(degree)),1));

e = int32(1);
nv = sfemesh_nnodes(mesh);
nhex = sfemesh_nelems(mesh,e);

mesh_buff_ = sfemesh_create(int32(3),nhex);
mesh_buff_.coords = mesh.coords;
keep_ = false(nhex,1);

% finding all hexes who share a face with a pyramid and keeping them
nkeep = int32(0);
for ii = 1:nhex
    [~, geid] = sfemesh_leid2teid(e, ii, mesh);
    jj = 1;
    while jj <= 6 && ~keep_(ii)
        hfid = mesh.sibhfs(geid,jj);
        if hfid
            fid = sfemesh_hfid2eid(hfid);
            [etable, ~] = sfemesh_teid2leid(mesh, fid);
            etype = mesh.elemtables(etable).etype;
            shape = obtain_elemshape(etype);
            if shape == SFE_SHAPE_PYRA
                keep_(ii) = true;
                nkeep = nkeep+1;
            end
        end
        jj = jj+1;
    end
end

% for every next layer find hexes that share a face with any hex in keep
if nlayersbig > 1
    keep2_ = false(nhex,1);
    for L = 2:nlayersbig
        for ii = 1:nhex
            [~, geid] = sfemesh_leid2teid(e, ii, mesh);
            jj = 1;
            while jj <= 6 && ~keep2_(ii) && ~keep_(ii)
                hfid = mesh.sibhfs(geid,jj);
                if hfid
                    fid = sfemesh_hfid2eid(hfid);
                    [etable, eid] = sfemesh_teid2leid(mesh, fid);
                    etype = mesh.elemtables(etable).etype;
                    shape = obtain_elemshape(etype);
                    if shape == SFE_SHAPE_HEXA && keep_(eid)
                        keep2_(ii) = true;
                        nkeep = nkeep+1;
                    end
                end
                jj=jj+1;
            end
        end

        for ii = 1:nhex
            keep_(ii) = keep_(ii) || keep2_(ii);
        end
    end
end

% adding kept hexahedrals to the buffer mesh
nelems = sfemesh_nelems(mesh,e);
etype = mesh.elemtables(e).etype;
nnodes = obtain_elemnnodes(etype);
[mesh_buff_,et] = sfemesh_append_econntable(mesh_buff_, etype);
mesh_buff_ = sfemesh_resize_econn(mesh_buff_,et,nkeep);
nn = int32(0);
for ii = 1:nelems
    if keep_(ii)
        nn=nn+1;
        for jj = 1:nnodes
            mesh_buff_.elemtables(et).conn(nn,jj) = mesh.elemtables(e).conn(ii,jj);
        end
    end
end


% adding tets and pyamids to buffer mesh
for etable = 2:sfemesh_nelemtypes(mesh)
    nelems = sfemesh_nelems(mesh,etable);
    etype = mesh.elemtables(etable).etype;
    nnodes = obtain_elemnnodes(etype);
    [mesh_buff_,et] = sfemesh_append_econntable(mesh_buff_, etype);
    mesh_buff_ = sfemesh_resize_econn(mesh_buff_,et,nelems);
    for ii = 1:nelems
        for jj = 1:nnodes
            mesh_buff_.elemtables(et).conn(ii,jj) = mesh.elemtables(etable).conn(ii,jj);
        end
    end
end

% decomposing the mesh to a linear one
mesh_buff_ = sfemesh_decompose_hiorder3d(mesh_buff_);

% finding hex layers in linear mesh
nhex = sfemesh_nelems(mesh_buff_,1);
keep_ = false(nhex,1);
nkeep = int32(0);
e = int32(1);
for ii = 1:nhex
    [~, geid] = sfemesh_leid2teid(e, ii, mesh_buff_);
    jj = 1;
    while jj <= 6 && ~keep_(ii)
        hfid = mesh_buff_.sibhfs(geid,jj);
        if hfid
            fid = sfemesh_hfid2eid(hfid);
            [etable, ~] = sfemesh_teid2leid(mesh_buff_, fid);
            etype = mesh_buff_.elemtables(etable).etype;
            shape = obtain_elemshape(etype);
            if shape == SFE_SHAPE_PYRA
                keep_(ii) = true;
                nkeep = nkeep+1;
            end
        end
        jj = jj+1;
    end
end

% for every next layer find hexes that share a face with any hex in keep
if nlayers > 1
    keep2_ = false(nhex,1);
    for L = 2:nlayers
        for ii = 1:nhex
            [~, geid] = sfemesh_leid2teid(e, ii, mesh_buff_);
            jj = 1;
            while jj <= 6 && ~keep2_(ii) && ~keep_(ii)
                hfid = mesh_buff_.sibhfs(geid,jj);
                if hfid
                    fid = sfemesh_hfid2eid(hfid);
                    [etable, eid] = sfemesh_teid2leid(mesh_buff_, fid);
                    etype = mesh_buff_.elemtables(etable).etype;
                    shape = obtain_elemshape(etype);
                    if shape == SFE_SHAPE_HEXA && keep_(eid)
                        keep2_(ii) = true;
                        nkeep = nkeep+1;
                    end
                end
                jj=jj+1;
            end
        end

        for ii = 1:nhex
            keep_(ii) = keep_(ii) || keep2_(ii);
        end
    end
end


% removing hexes that arent kept
nn = int32(0);
for ii = 1:nhex
    if keep_(ii)
        nn=nn+1;
        for jj = 1:8
            mesh_buff_.elemtables(e).conn(nn,jj) = mesh_buff_.elemtables(e).conn(ii,jj);
        end
    end
end
mesh_buff_ = sfemesh_resize_econn(mesh_buff_,e,nkeep);

% removing unused nodes
[mesh_buff_,mesh2aes,aes2mesh] = remove_unused_nodes(mesh_buff_);
nv = sfemesh_nnodes(mesh_buff_);

if coder.target("MATLAB")
    writevtk_unstr('mlin_1.vtk',mesh_buff_.coords, mesh_buff_.elemtables(1).conn);
    writevtk_unstr('mlin_2.vtk',mesh_buff_.coords, mesh_buff_.elemtables(2).conn);
    writevtk_unstr('mlin_3.vtk',mesh_buff_.coords, mesh_buff_.elemtables(3).conn);
    writevtk_unstr('mlin_4.vtk',mesh_buff_.coords, mesh_buff_.elemtables(4).conn);
end

% creating aes_mesh
aes_mesh = sfemesh_create(int32(3),nv);
aes_mesh.coords = mesh_buff_.coords;

% finding nelems from buffer mesh
nelems = int32(0);
ne = sfemesh_nelemtypes(mesh_buff_);
for e = 1:ne
    shape = obtain_elemshape(mesh_buff_.elemtables(e).etype);
    if shape == SFE_SHAPE_TET
        nelems = nelems + sfemesh_nelems(mesh_buff_,e);
    elseif shape == SFE_SHAPE_PYRA
        nelems = nelems + 2*sfemesh_nelems(mesh_buff_,e);
    elseif shape == SFE_SHAPE_HEXA
        nelems = nelems + 6*sfemesh_nelems(mesh_buff_,e);
    end
end

aes_mesh = sfemesh_append_econntable(aes_mesh,SFE_TET_4);
aes_mesh = sfemesh_resize_econn(aes_mesh,int32(1),nelems);
hex2tets = int32([1, 2, 4, 6; 2, 3, 4, 6; 5, 8, 6, 4; 1, 6, 4, 5; 3, 4, 6, 7; 4, 7, 8, 6]);
pyr2tets = int32([1,2,3,5, 1,3,4,5; 1,2,4,5, 2,3,4,5]);

hexfacesplits = int32([2,1,2,1; 1,2,1,2; 2,1,2,1; 2,1,2,1; 2,1,2,1; 2,1,2,1]);

eid = int32(0);
for e = 1:ne
    nelems = sfemesh_nelems(mesh_buff_,e);
    shape = obtain_elemshape(mesh_buff_.elemtables(e).etype);
    if shape == SFE_SHAPE_TET
        for ii = 1:nelems
            eid = eid+1;
            for jj = int32(1):4
                aes_mesh.elemtables(1).conn(eid,jj) = mesh_buff_.elemtables(e).conn(ii,jj);
            end
        end
    elseif shape == SFE_SHAPE_PYRA
        done_ = false(nelems,1);
        for ii = 1:nelems
            if ~done_(ii)
                [~, geid] = sfemesh_leid2teid(e, ii, mesh_buff_);
                hfid = mesh_buff_.sibhfs(geid,1);
                opp = sfemesh_hfid2eid(hfid);
                [e2, opp_eid] = sfemesh_teid2leid(mesh_buff_, opp);
                shape2 = obtain_elemshape(mesh_buff_.elemtables(e2).etype);
                lid = sfemesh_hfid2lid(hfid);
                if shape2 == SFE_SHAPE_HEXA
                    [~,leids] = obtain_facets(SFE_HEXA_8, int8(lid));
                    n = int32(0);
                    for jj = int32(1):4
                        if mesh_buff_.elemtables(e).conn(ii,1) == mesh_buff_.elemtables(e2).conn(opp_eid,leids(jj))
                            n = jj;
                        end
                    end

                    split = hexfacesplits(lid,n);
                    for kk = int32(1):2
                        eid = eid+1;
                        for jj = int32(1):4
                            aes_mesh.elemtables(1).conn(eid,jj) = mesh_buff_.elemtables(e).conn(ii, pyr2tets(split,4*(kk-1)+jj));
                        end
                    end
                    done_(ii) = true;
                else
                    m2cAssert(shape2 == SFE_SHAPE_PYRA);
                    [~,leids] = obtain_facets(SFE_HEXA_8, int8(lid));
                    n = int32(0);
                    for jj = int32(1):4
                        if mesh_buff_.elemtables(e).conn(ii,1) == mesh_buff_.elemtables(e2).conn(opp_eid,leids(jj))
                            n = jj;
                        end
                    end

                    for kk = int32(1):2
                        eid = eid+1;
                        for jj = int32(1):4
                            aes_mesh.elemtables(1).conn(eid,jj) = mesh_buff_.elemtables(e).conn(ii, pyr2tets(1,4*(kk-1)+jj));
                        end
                    end
                    done_(ii) = true;

                    for kk = int32(1):2
                        eid = eid+1;
                        for jj = int32(1):4
                            aes_mesh.elemtables(1).conn(eid,jj) = mesh_buff_.elemtables(e2).conn(opp_eid, pyr2tets((1-mod(n,2))+1,4*(kk-1)+jj));
                        end
                    end
                    done_(opp_eid) = true;
                end
            end
        end
    elseif shape == SFE_SHAPE_HEXA
        for ii = 1:nelems
            for kk = int32(1):6
                eid = eid+1;
                for jj = int32(1):4
                    aes_mesh.elemtables(1).conn(eid,jj) = mesh_buff_.elemtables(e).conn(ii, hex2tets(kk,jj));
                end
            end
        end
    end
end

aes_mesh = sfemesh_setup(aes_mesh,int32(-1),true);

% this part plots boundary tris to check the mesh is correct
%{
nelems = sfemesh_nelems(aes_mesh,1);
tris = zeros(nelems,3,'int32');
nt = int32(0);
for ii = 1:nelems
    for jj = 1:4
        if ~aes_mesh.sibhfs(ii,jj)
            [~,leids] = obtain_facets(SFE_TET_4, int8(jj));
            nt = nt+1;
            tris(nt,1:3) = aes_mesh.elemtables(1).conn(ii,leids'); 
        end
    end
end
tris = tris(1:nt,1:3);
writevtk_unstr('aes_t.vtk',aes_mesh.coords, tris);
%}
if coder.target("MATLAB")
    writevtk_unstr('aes.vtk',aes_mesh.coords, aes_mesh.elemtables(1).conn);
end

end