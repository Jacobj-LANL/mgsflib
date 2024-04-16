function [mesh,map,map_rev] = remove_unused_nodes(mesh)

%#codegen -args {SfeMesh}

coder.inline('never');

nv = int32(size(mesh.coords,1));
ndim = int32(size(mesh.coords,2));
coords_buff_ = zeros(nv,ndim);
map = zeros(nv,1,'int32');
map_rev = zeros(nv,1,'int32');
netypes = sfemesh_nelemtypes(mesh);
n = int32(0);
for e = int32(1):netypes
    nelems = sfemesh_nelems(mesh, e);
    nnodes = obtain_elemnnodes(mesh.elemtables(e).etype);
    elems_buff_ = zeros(nelems,nnodes,'int32');
    for ii = 1:nelems
        for jj = 1:nnodes
            v = mesh.elemtables(e).conn(ii,jj);
            if ~map(v)
                n = n+1;
                map(v) = n;
                map_rev(n) = v;
                coords_buff_(n,1:ndim) = mesh.coords(v,1:ndim);
            end
            elems_buff_(ii,jj) = map(v);
        end
    end

    mesh.elemtables(e).conn = elems_buff_;
end

mesh = sfemesh_resize_coords(mesh,n);
mesh.coords = coords_buff_(1:n,1:ndim);
map_rev = map_rev(1:n);
mesh = sfemesh_setup(mesh,int32(-1),true);
end