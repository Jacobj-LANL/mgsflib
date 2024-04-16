function [mesh,aes_mesh,mesh2aes,aes2mesh] = mgsfmesh_aesmesh_setup_tables(mesh,aes_mesh,mesh2aes,aes2mesh,nnsets)

%#codegen -args {SfeMesh, SfeMesh, coder.typeof(int32(0), [inf, 1]), coder.typeof(int32(0), [inf, 1]), coder.typeof(int32(0), [inf, 1])}

if coder.target('MATLAB') && exist(['meshing_gateway.' mexext], 'file')
    [mesh,aes_mesh,mesh2aes,aes2mesh] = meshing_gateway('mgsfmesh_aesmesh_setup_tables', mesh,aes_mesh,mesh2aes,aes2mesh,nnsets);
    return
end

nsets = length(nnsets);
nv_aes = sfemesh_nnodes(aes_mesh);
fnodetags_ = false(nv_aes,nsets);
xs_normals = zeros(nv_aes,3);
for nn = 1:nsets
    nsid = nnsets(nn);
    for ii = 1:sfemesh_nnodesets(mesh,nsid)
        fnodetags_(mesh2aes(mesh.nodesets(nsid).nids(ii)),nn) = true;
        xs_normals(mesh2aes(mesh.nodesets(nsid).nids(ii)),1:3) = mesh.nrmstables(nsid).normals(ii,1:3);
    end
end

[aes_mesh,fset] = sfemesh_append_facetset(aes_mesh,'f');
aes_mesh = sfemesh_resize_facetsetdata(aes_mesh,fset,nv_aes);
for nn = int32(1):nsets
   [aes_mesh,fset] = sfemesh_append_facetset(aes_mesh,sprintf('f%d',nn));
   aes_mesh = sfemesh_resize_facetsetdata(aes_mesh,fset,nv_aes);
end
nfs = zeros(1,nsets,'int32');
nf1 = int32(0);

faces = int32([1,3,2; 1,2,4; 2,3,4; 3,1,4]);
for ii = 1:sfemesh_nelems(aes_mesh)
    for jj = int32(1):4
        if ~aes_mesh.sibhfs(ii,jj)
            v1 = aes_mesh.elemtables(1).conn(ii,faces(jj,1));
            v2 = aes_mesh.elemtables(1).conn(ii,faces(jj,2));
            v3 = aes_mesh.elemtables(1).conn(ii,faces(jj,3));
            found = false;
            for nn = 1:nsets
                if fnodetags_(v1,nn) && fnodetags_(v2,nn) && fnodetags_(v3,nn)
                    nfs(nn) = nfs(nn)+1;
                    aes_mesh.facetsets(nn+1).eids(nfs(nn)) = ii;
                    aes_mesh.facetsets(nn+1).facets(nfs(nn)) = jj;
                    found = true;
                end
            end

            if ~found
                nf1 = nf1+1;
                aes_mesh.facetsets(1).eids(nf1) = ii;
                aes_mesh.facetsets(1).facets(nf1) = jj;
            end
        end
    end
end

aes_mesh = sfemesh_resize_facetsetdata(aes_mesh,int32(1),nf1);
for nn = 1:nsets
   aes_mesh = sfemesh_resize_facetsetdata(aes_mesh,int32(nn+1),nfs(nn));
end

s = 'f';
for nn = int32(1):nsets
    s2 = sprintf(',f%d',nn);
    s = [s,s2];
end
aes_mesh = mgsfmesh_facetsets2nodesets(aes_mesh, s, ...
        int32([100; 200*ones(nsets,1)]), int32([100; 200]));

aes_mesh = sfemesh_nodefacetsets2nrmstables(aes_mesh, int32([100; 200]), xs_normals);
end