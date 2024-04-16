function mesh = create_sphere_meshmized_nrmstables(mesh)
normal_imp = @(x,y,z) -[2*x-1, 2*y-1, 2*z-1];

% creating normals arrays
xs_normals = zeros(sfemesh_nnodes(mesh), 3);
for ii = 1:length(mesh.nodesets(2).nids)
    xs_normals(mesh.nodesets(2).nids(ii), :) = normal_imp( ...
        mesh.coords(mesh.nodesets(2).nids(ii), 1), ...
        mesh.coords(mesh.nodesets(2).nids(ii), 2), ...
        mesh.coords(mesh.nodesets(2).nids(ii), 3));
end

% calculating remaining normals and normalizing
normals = [0, 0, -1; 0, -1, 0; 1, 0, 0; ...
    0, 1, 0; 0, 0, 1; 0, -1, 0];
for nf = int32(1):6
    for ii = 1:length(mesh.facetsets(nf).eids)
        geid = mesh.facetsets(nf).eids(ii);
        facet = mesh.facetsets(nf).facets(ii);
        [etable, eid] = sfemesh_teid2leid(mesh, geid);
        etype = mesh.elemtables(etable).etype;
        [~, leids] = obtain_facets(etype, facet);
        for jj = 1:length(leids)
            xs_normals(mesh.elemtables(etable).conn(eid, leids(jj)), :) = ...
                xs_normals(mesh.elemtables(etable).conn(eid, leids(jj)), :) + normals(nf, :);
        end
    end
end

for ii = 1:size(xs_normals, 1)
    if abs(norm(xs_normals(ii,:))) > 1e-6
        xs_normals(ii, :) = xs_normals(ii, :) ./ norm(xs_normals(ii, :));
    end
end

% creating normalset data
mesh = sfemesh_nodefacetsets2nrmstables(mesh, int32([100; 100; 100; 100; 100; 100; 200]), xs_normals);
end