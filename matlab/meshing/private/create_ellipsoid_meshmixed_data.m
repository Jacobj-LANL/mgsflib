function mesh = create_ellipsoid_meshmixed_data(mesh,n)
% setup mixed mesh for ellipsoidhole geometry
normal_imp = @(x,y,z) -[(8*x)/9 - 4/9, 2*y - 1, (32*z)/9 - 16/9];

tol = 1e-8;

% create Left, Right, Top, Bottom facet groups
mesh = sfemesh_append_facetset(mesh, 'Bottom');
mesh = sfemesh_resize_facetsetdata(mesh, int32(1), (n - 1) * (n - 1));
bottom = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Front');
mesh = sfemesh_resize_facetsetdata(mesh, int32(2), (n - 1) * (n - 1));
front = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Right');
mesh = sfemesh_resize_facetsetdata(mesh, int32(3), (n - 1) * (n - 1));
right = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Back');
mesh = sfemesh_resize_facetsetdata(mesh, int32(4), (n - 1) * (n - 1));
back = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Top');
mesh = sfemesh_resize_facetsetdata(mesh, int32(5), (n - 1) * (n - 1));
top = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Left');
mesh = sfemesh_resize_facetsetdata(mesh, int32(6), (n - 1) * (n - 1));
left = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Surface');
mesh = sfemesh_resize_facetsetdata(mesh, int32(7), sfemesh_nelems(mesh,1));
surface = int32(0);

% creating facetsets for the exterior box
nelems = sfemesh_nelems(mesh,2);
istart = mesh.elemtables(2).istart;
etype = mesh.elemtables(2).etype;
for e = 1:nelems
    eid = e + istart - 1;
    for j = 1:6
        if mesh.sibhfs(eid, j) ~= 0; continue; end
        [ftype, lids] = obtain_facets(etype, int8(j));
        xsloc_ = m2cNullcopy(zeros(4, 3));
        v_ = m2cNullcopy(zeros(4, 1, 'int32'));
        for lid = 1:4
            v_(lid) = mesh.elemtables(2).conn(e, lids(lid));
        end
        for lid = 1:4
            xsloc_(lid, 1) = mesh.coords(v_(lid), 1);
            xsloc_(lid, 2) = mesh.coords(v_(lid), 2);
            xsloc_(lid, 3) = mesh.coords(v_(lid), 3);
        end

        [N, ~] = normal3d(xsloc_);
        if abs(N(3) + 1) < tol
            bottom = bottom + 1;
            mesh.facetsets(1).eids(bottom) = eid;
            mesh.facetsets(1).facets(bottom) = j;

        elseif abs(N(1) - 1) < tol
            right = right + 1;
            mesh.facetsets(3).eids(right) = eid;
            mesh.facetsets(3).facets(right) = j;
        elseif abs(N(3) - 1) < tol
            top = top + 1;
            mesh.facetsets(5).eids(top) = eid;
            mesh.facetsets(5).facets(top) = j;

        elseif abs(N(1) + 1) < tol
            left = left + 1;
            mesh.facetsets(6).eids(left) = eid;
            mesh.facetsets(6).facets(left) = j;

        elseif abs(N(2) + 1) < tol
            front = front + 1;
            mesh.facetsets(2).eids(front) = eid;
            mesh.facetsets(2).facets(front) = j;

        elseif abs(N(2) - 1) < tol
            back = back + 1;
            mesh.facetsets(4).eids(back) = eid;
            mesh.facetsets(4).facets(back) = j;
        end
    end
end

% creating facetsets for interior surface(s)
nelems = sfemesh_nelems(mesh,1);
istart = mesh.elemtables(1).istart;
etype = mesh.elemtables(1).etype;
for e = 1:nelems
    eid = e + istart - 1;
    for j = 1:4
        if mesh.sibhfs(eid, j) ~= 0; continue; end
        surface = surface + 1;
        mesh.facetsets(7).eids(surface) = eid;
        mesh.facetsets(7).facets(surface) = j;
    end
end

% resizing
mesh = sfemesh_resize_facetsetdata(mesh, int32(1), bottom);
mesh = sfemesh_resize_facetsetdata(mesh, int32(2), front);
mesh = sfemesh_resize_facetsetdata(mesh, int32(3), right);
mesh = sfemesh_resize_facetsetdata(mesh, int32(4), back);
mesh = sfemesh_resize_facetsetdata(mesh, int32(5), top);
mesh = sfemesh_resize_facetsetdata(mesh, int32(6), left);
mesh = sfemesh_resize_facetsetdata(mesh, int32(7), surface);

% creating nodesets
mesh = sfemesh_facetsets2nodesets(mesh, 'Bottom,Front,Right,Back,Top,Left,Surface', ...
    [100; 100; 100; 100; 100; 100; 200], [100; 200]);

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