function [mesh,tetmesh] = mgsfmesh_shape_mixed(shape, n, degree, dist, hgr, regenerate, min_angle)
%mgsfmesh_shape_mixed - Generate mixed mgsfmesh with shape hole.
%
%   mesh = mgsfmesh_shape_mixed(shape, n, degree, dist, hgr, regenerate)
%   mesh = mgsfmesh_shape_mixed(shape, n, degree, dist, hgr)
%   mesh = mgsfmesh_shape_mixed(shape, n, degree)
%   mesh = mgsfmesh_shape_mixed(shape, n)
%   mesh = mgsfmesh_shape_mixed(shape)
%
% Parameters
% -----------
%   shape: name of the .geo model to be used (must be in meshing/private/models)
%   n: number of points in x,y,z-direction (int32)
%   degree: degree of mesh
%   dist: nodal distribution of high order nodes
%   hgr: ratio of minimum size of boundary elems to interior elems
%   regenerate: whether or not to regenerate mesh or used saved linear mesh
%
% Returns
% --------
%   mesh: an SfeMesh object.
%
% Notes
% -----
% mesh generated using Gmsh and tetgen

if nargin < 7
    min_angle = 10;
end
if nargin < 6
    regenerate = true;
end
if nargin < 5
    hgr = 0.1;
end
if nargin < 4
    dist = 'Gauss-Lobatto';
end
if nargin < 3
    degree = 3;
end
if nargin < 2
    n = 25;
end
if nargin < 1
    shape = 'sphere';
end

% if linear mesh exists load it
filename = [mgsflib_root,'/meshing/saved_meshes/',shape, '_',num2str(n),'_',num2str(hgr~=1),'_linear.mat'];
if ~regenerate && isfile(filename)
    load(filename);
else

tic;
% Use gmsh to create the h-gr surface mesh
modelname = [mgsflib_root,'/meshing/private/models/',shape];
[surf_tris, surf_xs, surf_normals, inflags] = gmsh_surface(string(modelname), int32(n), hgr);
for ii = 1:size(surf_normals,1)
    surf_normals(ii,2:4) = surf_normals(ii,2:4)/norm(surf_normals(ii,2:4));
end

%{
figure;
hold on
%axis([0.3,0.7,0.3,0.7,0.3,0.7])
trimesh(surf_tris,surf_xs(:,1),surf_xs(:,2),surf_xs(:,3),'EdgeColor','k')
quiver3(surf_xs(surf_normals(:,1),1),surf_xs(surf_normals(:,1),2),surf_xs(surf_normals(:,1),3),surf_normals(:,2),surf_normals(:,3),surf_normals(:,4))
axis equal
%}

% if excplicit formula is given for determining inside / outside, use it
if exist('explicit')
    bar = linspace(0, 1, n);
    inflags = false(n * n * n, 1);
    p = int32(0);
    for k = 1:n
        for j = 1:n
            for i = 1:n
                p = p + 1;
                if explicit(bar(i), bar(j), bar(k)) < 0
                    inflags(p) = true;
                end
            end
        end
    end
end

% creating staircase mesh and extracting boundary
mesh = sfemesh_create_staircase_mesh(inflags, int32([n,n,n]), [0,0,0;1,1,1]);
nv = sfemesh_nnodes(mesh);
bnd_mesh = extract_bndmesh3d(mesh, 1);

% merging two surface meshes for tetgen input
nvbnd = sfemesh_nnodes(bnd_mesh);
nebnd = sfemesh_nelems(bnd_mesh);
nvtotal = nvbnd + size(surf_xs,1);
netotal = nebnd + size(surf_tris,1);
bnd_mesh = sfemesh_resize_coords(bnd_mesh, nvtotal);
bnd_mesh.coords(nvbnd + 1:nvtotal, 1:3) = surf_xs;
bnd_mesh = sfemesh_resize_econn(bnd_mesh, 1, netotal);
bnd_mesh.elemtables(1).conn(nebnd + 1:netotal, 1:3) = surf_tris + nvbnd;

% generating interior hole using tetgen
inputstring = "-pqk1.2AaY";
h = 1/(double(n)-1);
targetvol = 1/6 * h^3;
[tets, tet_coords] = tetgen(inputstring, bnd_mesh.elemtables(1).conn, bnd_mesh.coords,[0.5, 0.5, 0.5 + 0.2+h/4, targetvol], [0.5,0.5,0.5]);

% merging tets into mesh
nvtet = size(tet_coords,1);
ntets_hole = size(tets,1);
nv_new = nv + nvtet - nvbnd;
mesh = sfemesh_resize_coords(mesh, nv_new);
mesh.coords((nv + 1):nv_new, 1:3) = tet_coords(nvbnd + 1:nvtet, 1:3);
nelems_mixed = sfemesh_nelems(mesh, 3);
mesh = sfemesh_resize_econn(mesh, 3, ntets_hole + nelems_mixed);
nbnd = int32(length(bnd_mesh.nodesets(1).nids));
for ii = 1:ntets_hole
    for jj = int32(1):4
        if tets(ii, jj) > nbnd
            mesh.elemtables(3).conn(ii + nelems_mixed, jj) = tets(ii, jj) + nv - nvbnd;
        else
            mesh.elemtables(3).conn(ii + nelems_mixed, jj) = bnd_mesh.nodesets(1).nids(tets(ii, jj));
        end
    end
end

% call stellar for tet only part
mesh = mgsfmesh_improve_tets_stellar_mixed(mesh);
mesh.facetsets = m2cNullcopy(repmat(FacetSet(char(zeros(1, m2cZero))), m2cZero, 1));
disp(['finished generating linear mesh in ', num2str(toc),' seconds']);

save(filename,"mesh");

end

if min_angle ~= 0
    valid_eids = zeros(100,1,'int32');
    neids = int32(0);
    for ii = 1:size(mesh.elemtables(3).conn,1)
        eid = mesh.elemtables(3).istart + ii - 1;
        jj = 1; stop = false;
        while jj <= 4 && ~stop
            hfid = mesh.sibhfs(eid,jj);
            if hfid
                if sfemesh_hfid2eid(hfid) >= mesh.elemtables(3).istart
                    neids = neids + 1;
                    valid_eids(neids) = ii;
                    stop = true;
                end
            end
            jj=jj+1;
        end
    end
    eids = valid_eids(randi(neids));
    eids(1) = 2018;
    
    
    angles_elem = zeros(1,6);
    for ii = 1:length(eids)
        eid = eids(ii);
        normal = mesh.coords(mesh.elemtables(3).conn(eid,4),1:3)-mesh.coords(mesh.elemtables(3).conn(eid,3),1:3);

        xs = mesh.coords(mesh.elemtables(3).conn(eid,1:4),1:3);
        P1 = xs(1,:);
        P2 = xs(2,:);
        P3 = xs(3,:);
        P4 = xs(4,:);
   
        angles_elem(1) = dihedral(P2-P1,P3-P1,P4-P1);
        angles_elem(2) = dihedral(P3-P2,P1-P2,P4-P2);
        angles_elem(3) = dihedral(P1-P3,P2-P3,P4-P3);
        angles_elem(4) = dihedral(P4-P1,P2-P1,P3-P1);
        angles_elem(5) = dihedral(P4-P2,P1-P2,P3-P2);
        angles_elem(6) = dihedral(P4-P3,P1-P3,P2-P3); 
        
        mesh.coords(mesh.elemtables(3).conn(eid,4),1:3) = mesh.coords(mesh.elemtables(3).conn(eid,3),1:3) ...
            + (min_angle/angles_elem(1))*normal/1.27303;
    end
end

posjac = sfemesh_check_jacobians(mesh);
angles = min_dihedral_angles(mesh);
fprintf(1,'minimum dihedral angle in linear mesh: %g\n',min(angles));


if false
    writevtk_unstr('mesh_lin1.vtk',mesh.coords, mesh.elemtables(1).conn);
    writevtk_unstr('mesh_lin2.vtk',mesh.coords, mesh.elemtables(2).conn);
    writevtk_unstr('mesh_lin3.vtk',mesh.coords, mesh.elemtables(3).conn);
end

if strcmp('Gauss-Lobatto',dist)
    filename_high = [mgsflib_root,'/meshing/saved_meshes/',shape, '_',num2str(n),'_',num2str(hgr~=1),'_',num2str(degree),'_gl.mat'];
else
    filename_high = [mgsflib_root,'/meshing/saved_meshes/',shape, '_',num2str(n),'_',num2str(hgr~=1),'_',num2str(degree),'_eq.mat'];
end
if ~regenerate && isfile(filename_high)
    load(filename_high); 
    posjac = sfemesh_check_jacobians(mesh);
    if ~posjac
        dbstop;
    end
    return;
end


if degree > 1
    mesh = sfemesh_raise_hiorder3d(mesh, int32(degree), dist);
    if strcmp(shape,'sphere')
        mesh = create_sphere_meshmixed_data(mesh,n);
        [mesh, tetmesh] = mgsfmesh_extract_tets(mesh, int32(min(degree+2, 6)), int32(degree+2 < 5 && obtain_elemnodepos(mesh.elemtables(3).etype)),int32(2), int32(7));
        [mesh, tetmesh] = mgsfmesh_raise_curved_sphere(mesh,tetmesh,0.05);
        mesh = create_sphere_meshmized_nrmstables(mesh);
        posjac = sfemesh_check_jacobians(tetmesh);
        save(filename_high,"mesh","tetmesh");
    else
        save(filename_high,"mesh");
    end
    posjac = sfemesh_check_jacobians(mesh);
else
    if strcmp(shape,'sphere')
        mesh = create_sphere_meshmixed_data(mesh,n);
        save(filename_high,"mesh");
    end
end


 
    %{
    figure;
    hold on
    for ii = 1:sfemesh_nfacetsets(mesh,int32(7))
        clf;
        geid = mesh.facetsets(7).eids(ii);
        [~, eid] = sfemesh_teid2leid(mesh, geid);
        xyz = tetmesh.coords(tetmesh.elemtables(1).conn(eid,:),1:3);
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'filled');
        pause(0.1);
    end

    for ii = 1:sfemesh_nedgesets(mesh,int32(1))
        clf;
        eid = mesh.edgesets(1).eids(ii);
        xyz = tetmesh.coords(tetmesh.elemtables(1).conn(eid,:),1:3);
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'filled');
        pause(1);
    end
    %}

end

function angles = min_dihedral_angles(mesh)
    nelems = sfemesh_nelems(mesh,3);
    angles = zeros(nelems,1);
    angles_elem = zeros(1,6);
    for ii = 1:nelems
        xs = mesh.coords(mesh.elemtables(3).conn(ii,1:4),1:3);
        P1 = xs(1,:);
        P2 = xs(2,:);
        P3 = xs(3,:);
        P4 = xs(4,:);
   
        angles_elem(1) = dihedral(P2-P1,P3-P1,P4-P1);
        angles_elem(2) = dihedral(P3-P2,P1-P2,P4-P2);
        angles_elem(3) = dihedral(P1-P3,P2-P3,P4-P3);
        angles_elem(4) = dihedral(P4-P1,P2-P1,P3-P1);
        angles_elem(5) = dihedral(P4-P2,P1-P2,P3-P2);
        angles_elem(6) = dihedral(P4-P3,P1-P3,P2-P3);  
        angles(ii) = min(angles_elem);
    end
end

function angle = dihedral(b0, b1, b2)
    u = cross(b0,b1);
    v = cross(b0,b2);
    
    angle = acos((u*v')/(norm(u)*norm(v)))*180/pi;
end

function [N] = normal3d(xsloc)

coder.inline('always')

    v1 = [xsloc(2, 1) - xsloc(1, 1), xsloc(2, 2) - xsloc(1, 2), xsloc(2, 3) - xsloc(1, 3)];
    v2 = [xsloc(3, 1) - xsloc(2, 1), xsloc(3, 2) - xsloc(2, 2), xsloc(3, 3) - xsloc(2, 3)];
    v3 = [xsloc(1, 1) - xsloc(3, 1), xsloc(1, 2) - xsloc(3, 2), xsloc(1, 3) - xsloc(3, 3)];

    N = [v1(2) * -v3(3) - v1(3) * -v3(2), v1(3) * -v3(1) - v1(1) * -v3(3), v1(1) * -v3(2) - v1(2) * -v3(1)];
end