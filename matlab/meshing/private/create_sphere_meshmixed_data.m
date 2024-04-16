function mesh = create_sphere_meshmixed_data(mesh,n)
% setup mixed mesh for ellipsoidhole geometry

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
nelems = sfemesh_nelems(mesh,1);
istart = mesh.elemtables(1).istart;
etype = mesh.elemtables(1).etype;
for e = 1:nelems
    eid = e + istart - 1;
    for j = 1:6
        if mesh.sibhfs(eid, j) ~= 0; continue; end
        [ftype, lids] = obtain_facets(etype, int8(j));
        xsloc_ = m2cNullcopy(zeros(4, 3));
        v_ = m2cNullcopy(zeros(4, 1, 'int32'));
        for lid = 1:4
            v_(lid) = mesh.elemtables(1).conn(e, lids(lid));
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
nelems = sfemesh_nelems(mesh,3);
istart = mesh.elemtables(3).istart;
etype = mesh.elemtables(3).etype;
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
mesh = mgsfmesh_facetsets2nodesets(mesh, 'Bottom,Front,Right,Back,Top,Left,Surface', ...
    int32([100; 100; 100; 100; 100; 100; 200]), int32([100; 200]));

end