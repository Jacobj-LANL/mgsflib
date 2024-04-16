function mesh = setup_boundary3D_hole(mesh, nx, ny, nz)
%setup_boundary3D_hole - Generate nodeset, facetsets and normalsets for 3D simple mesh
%
%   mesh = setup_boundary3D_hole(mesh,nx,ny, h_ratio, fnormals)
%
% Parameters
% ----------
%   mesh: an SfeMesh object.
%   nx: number of points in x-direction.
%   ny: number of points in y-direction.
%
% Returns
% -------
%   mesh: an SfeMesh object.
% See also: ahmesh_gen_rectmesh, ahmesh_gen_trimesh
%
% Notes: Sets up the facetsets and normalsets for a mesh with exterior
% boundary on unit square and interior boundary along some curve facetsets
% 1:4 are for the 4 sides of the unit square and facetset 5 is for the
% interior curve.

tol = 1e-8;

% create Left, Right, Top, Bottom facet groups
mesh = sfemesh_append_facetset(mesh, 'Bottom');
mesh = sfemesh_resize_facetsetdata(mesh, int32(1), 2 * (nx - 1) * (ny - 1));
bottom = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Front');
mesh = sfemesh_resize_facetsetdata(mesh, int32(2), 2 * (nx - 1) * (nz - 1));
front = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Right');
mesh = sfemesh_resize_facetsetdata(mesh, int32(3), 2 * (ny - 1) * (nz - 1));
right = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Back');
mesh = sfemesh_resize_facetsetdata(mesh, int32(4), 2 * (nx - 1) * (nz - 1));
back = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Top');
mesh = sfemesh_resize_facetsetdata(mesh, int32(5), 2 * (nx - 1) * (ny - 1));
top = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Left');
mesh = sfemesh_resize_facetsetdata(mesh, int32(6), 2 * (ny - 1) * (nz - 1));
left = int32(0);

mesh = sfemesh_append_facetset(mesh, 'Surface');
mesh = sfemesh_resize_facetsetdata(mesh, int32(7), sfemesh_nnodes(mesh));
surface = int32(0);

% array to keep track of facets connected to each node
for etable = 1:sfemesh_nelemtypes(mesh)
    ne = sfemesh_nelems(mesh, etable);
    etype = sfemesh_etype(mesh, etable, int32(1));
    nf = obtain_facets(etype); % total number of facets
    istart = mesh.elemtables(etable).istart;
    assert(nf <= 6 && nf >= 1);
    for e = 1:ne
        eid = e + istart - 1;
        for j = 1:nf
            if mesh.sibhfs(eid, j) ~= 0; continue; end
            [ftype, lids] = obtain_facets(etype, int8(j));
            nlids = obtain_elemnnodes(ftype);
            xsloc_ = m2cNullcopy(zeros(nlids, 3));
            v_ = m2cNullcopy(zeros(nlids, 1, 'int32'));
            for lid = 1:nlids
                v_(lid) = mesh.elemtables(1).conn(eid, lids(lid));
            end
            for lid = 1:nlids
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

            else
                surface = surface + 1;
                mesh.facetsets(7).eids(surface) = eid;
                mesh.facetsets(7).facets(surface) = j;

            end
        end
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
end
