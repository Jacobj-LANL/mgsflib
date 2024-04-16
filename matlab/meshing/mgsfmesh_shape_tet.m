function [mesh,tetmesh] = mgsfmesh_shape_tet(shape, n, degree, dist, hgr, regenerate)
%mgsfmesh_shape_tet - Generate tet mgsfmesh with shape hole.
%
%   mesh = mgsfmesh_shape_tet(shape, n, degree, dist, hgr, regenerate)
%   mesh = mgsfmesh_shape_tet(shape, n, degree, dist, hgr)
%   mesh = mgsfmesh_shape_tet(shape, n, degree)
%   mesh = mgsfmesh_shape_tet(shape, n)
%   mesh = mgsfmesh_shape_tet(shape)
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
%   mesh: an mgsfmesh object.
%
% Notes
% -----
% mesh generated using Gmsh

if nargin < 6
    regenerate = true;
end
if nargin < 5
    hgr = 1;
end
if nargin < 4
    dist = 'Gauss-Lobatto';
end
if nargin < 3
    degree = 2;
end
if nargin < 2
    n = 20;
end
if nargin < 1
    shape = 'sphere';
end

% if linear mesh exists load it
filename = [mgsflib_root,'/meshing/saved_meshes/',shape, '_',num2str(n),'_',num2str(hgr~=1),'_tet_linear.mat'];
if ~regenerate && isfile(filename)
    load(filename);
else

% Use gmsh to create the h-gr surface mesh
modelname = [mgsflib_root,'/meshing/private/models/',shape];
tic;
[tets, xs] = gmsh_volume(string(modelname), n, hgr);

nv = size(xs,1);
mesh = sfemesh_create(int32(3),int32(nv));
mesh.coords = xs;

nelems = size(tets,1);
mesh = sfemesh_append_econntable(mesh, SFE_TET_4);
mesh = sfemesh_resize_econn(mesh,int32(1),nelems);
mesh.elemtables(1).conn = tets;
mesh = sfemesh_setup(mesh,int32(-1),true);

disp(['finished generating linear mesh in ', num2str(toc),' seconds']);

save(filename,"mesh");
end

posjac = sfemesh_check_jacobians(mesh);
angles = min_dihedral_angles(mesh);
fprintf(1,'minimum dihedral angle in linear mesh: %g\n',min(angles));

if strcmp('Gauss-Lobatto',dist)
    filename_high = [mgsflib_root,'/meshing/saved_meshes/',shape, '_',num2str(n),'_',num2str(hgr~=1),'_',num2str(degree),'_gl_tet.mat'];
else
    filename_high = [mgsflib_root,'/meshing/saved_meshes/',shape, '_',num2str(n),'_',num2str(hgr~=1),'_',num2str(degree),'_eq_tet.mat'];
end

if ~regenerate && isfile(filename_high)
    load(filename_high); 
    posjac = sfemesh_check_jacobians(mesh);
    if ~posjac
        %dbstop;
    end
    return;
end

if degree > 1
    if strcmp(shape,'sphere')
        mesh = setup_boundary3D_hole(mesh, n, n, n);
        mesh = sfemesh_raise_hiorder3d(mesh, int32(degree), dist);
        % creating nodesets
        mesh = mgsfmesh_facetsets2nodesets(mesh, 'Bottom,Front,Right,Back,Top,Left,Surface', ...
        int32([100; 100; 100; 100; 100; 100; 200]), int32([100; 200]));
        [mesh, tetmesh] = mgsfmesh_extract_tets(mesh, int32(5), int32(0),int32(2), int32(7));
        [mesh, tetmesh] = mgsfmesh_raise_curved_sphere(mesh,tetmesh,0.05);
        mesh = create_sphere_meshmized_nrmstables(mesh);
    end
    posjac = sfemesh_check_jacobians(mesh);
else
    mesh = create_sphere_meshmixed_data(mesh,n);
end

save(filename_high,"mesh","tetmesh");

end

function angles = min_dihedral_angles(mesh)
    nelems = sfemesh_nelems(mesh,1);
    angles = zeros(nelems,1);
    angles_elem = zeros(1,6);
    for ii = 1:nelems
        xs = mesh.coords(mesh.elemtables(1).conn(ii,1:4),1:3);
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