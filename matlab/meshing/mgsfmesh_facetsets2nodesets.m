function [mesh, nsetids] = mgsfmesh_facetsets2nodesets(mesh, fsetnames, ...
    bctags, bcpriority)
%sfemesh_facetsets2nodesets - Convert facet sets to node sets
%
%   [mesh, nsetids] = sfemesh_facetsets2nodesets(mesh, fsetnames, bctags, bcpriority)
%
% Parameters
% ----------
%   mesh:       An SfeMesh instance
%   fsetnames:  Facet set names
%   bctags:     Boundary condition tags wrt `fsetnames`
%   bcpriority: Boundary condition priority order
%
% Returns
% -------
%   mesh:       `mesh` must be passed by reference
%   nsetids:    Node set IDs in `mesh.nodesets` corresponding to newly
%               added sets for each BC type
%
% Notes
% -----
% Given a list of facet sets by their names and corresponding user-defined
% boundary condition tags, this function combines and converts facets sets
% into node sets wrt the total number of user-defined boundary conditions.
% In addition, the boundary condition priority defines the node priority in
% resulting node sets; for instance, if BC1 takes higher priority than BC2,
% and the corresponding facet sets of BC1 and BC2 are adjacent, the the
% nodes along the joined regions belongs to the node set corresponding to BC1.
%
% The boundary condition tags can take arbitrary integer values defined by
% the user, but they must be unique wrt BC types. Also, one tag cannot have
% multiple priority orders in `bcpriority`.
%
% Regarding `fsetnames`, the input should be a char string of all facet set
% names with delimiter comma (','). Moreover, spaces are not allowed before
% and after delimiters.
%
% Examples
% --------
% Assume we have four facet sets, namely 'Bottom', 'Top', 'Right', and
% 'Left'. They have corresponding boundary condition tags 100, 300, 200,
% and 200 thus total three boundary condition types. Therefore, we have
%
%   fsetnames = 'Bottom,Top,Right,Left';
%   bctags = int32([100;300;200;200]);
%
% In the output node sets, we want BC 100 and 300 have the highest and
% lowest priority, respectively. Hence,
%
%   bcpriority = int32([100;200;300]);
%
% For the resulting three node sets, BC100 have all nodes in facet
% 'Bottom', and BC200 have all nodes in 'Right' and 'Left' except for those
% in 'Bottom'. The remaining nodes in the four facet sets are in BC300.
%
% See also sfemesh_append_nodeset, sfemesh_append_facetset

%#codegen -args {SfeMesh, coder.typeof(char(0), [1 inf]),
%#codegen        coder.typeof(int32(0), [inf 1]),
%#codegen        coder.typeof(int32(0), [inf 1])}

if coder.target('MATLAB') && exist(['meshing_gateway.' mexext], 'file')
    mesh = meshing_gateway('mgsfmesh_facetsets2nodesets', mesh, fsetnames,bctags, bcpriority);
    return
end

coder.inline('never');

if isempty(mesh.teids)
    error('sfemesh_facetsets2nodesets:missingSetup', 'call setup on this mesh');
end

% Step I: extract out corresponding facet set indices
nfacets = cast(size(bctags, 1), 'int32');
fsetids_ = facetsetnames2ids(mesh, fsetnames, nfacets);

% Step II: build priority list wrt facet indices in CRS format
[bc_ptr_, bc_fsets_] = build_bcpriority(bctags, bcpriority, fsetids_);

% Step III: convert facet sets to node sets wrt bc types
[mesh, nsetids] = build_nsets(mesh, bc_ptr_, bc_fsets_, bcpriority);

end

function fsetids = facetsetnames2ids(mesh, fsetnames, nfacets)
% convert facet set names to IDs

coder.inline('never');
m2cNowarnBuf('\w+'); % disable char buffer for error printing

if isempty(fsetnames)
    error('facetsetnames2ids:emptySetName', 'facet name set cannot be empty');
end

% determine how many names in fsetnames
nchar = cast(size(fsetnames, 2), 'int32');
i = int32(1);
nfacets1 = int32(0);
while i <= nchar
    if fsetnames(i) == ','
        nfacets1 = nfacets1 + 1;
        i = i + 1;
    end
    i = i + 1;
end
% add one if the last char is not deliminter
nfacets1 = nfacets1 + int32(fsetnames(end) ~= ',');

if nfacets1 ~= nfacets
    error('facetsetnames2ids:unmatchedSetSize', ...
        'facet size (%d) does not match input %d', nfacets1, nfacets);
end

fsetids = m2cNullcopy(zeros(nfacets1, 1, 'int32'));

% construct a map
map = m2cMapCreate;

i = int32(1);
p = int32(1);
while i <= nchar
    if fsetnames(i) == ','
        fsetname_ = m2cNullcopy(fsetnames(p:i - 1));
        for j = p:i - 1; fsetname_(j - p + 1) = fsetnames(j); end
        if m2cMapIsKey(map, fsetname_)
            map = m2cMapDestroy(map);
            error('facetsetnames2ids:dupNames', '%s facet set is duplicated', fsetname_);
        end
        map = m2cMapSetValue(map, fsetname_, int32(0));
        i = i + 1;
        p = i;
    end
    i = i + 1;
end
if fsetnames(end) ~= ','
    fsetname_ = m2cNullcopy(fsetnames(p:end));
    for j = p:nchar; fsetname_(j - p + 1) = fsetnames(j); end
    if m2cMapIsKey(map, fsetname_)
        map = m2cMapDestroy(map);
        error('facetsetnames2ids:dupNames', '%s facet set is duplicated', fsetname_);
    end
    map = m2cMapSetValue(map, fsetname_, int32(0));
end

% register facet set to the map
nfacetmesh = sfemesh_nfacetsets(mesh);
for i = 1:nfacetmesh
    if ~m2cMapIsKey(map, mesh.facetsets(i).name); continue; end
    map = m2cMapSetValue(map, mesh.facetsets(i).name, i);
end

i = int32(1);
p = int32(1);
fscount = int32(1);
while i <= nchar
    if fsetnames(i) == ','
        fsetname_ = m2cNullcopy(fsetnames(p:i - 1));
        for j = p:i - 1; fsetname_(j - p + 1) = fsetnames(j); end
        fsetid = m2cMapGetValue(map, fsetname_);
        if fsetid == 0
            map = m2cMapDestroy(map);
            error('facetsetnames2ids:badFasetSet', '%s facet not exist', fsetname_);
        end
        i = i + 1;
        p = i;
        fsetids(fscount) = fsetid;
        fscount = fscount + 1;
    end
    i = i + 1;
end
if fsetnames(end) ~= ','
    fsetname_ = m2cNullcopy(fsetnames(p:end));
    for j = p:nchar; fsetname_(j - p + 1) = fsetnames(j); end
    fsetid = m2cMapGetValue(map, fsetname_);
    if fsetid == 0
        map = m2cMapDestroy(map);
        error('facetsetnames2ids:badFasetSet', '%s facet not exist', fsetname_);
    end
    fsetids(end) = fsetid;
end

map = m2cMapDestroy(map); %#ok<*NASGU>

end

function [bc_ptr, bc_fsets] = build_bcpriority(bctags, bcpriority, fsetids)
% build priority list

coder.inline('never');

nfacets = cast(size(bctags, 1), 'int32');
nbcs = cast(size(bcpriority, 1), 'int32');

map = m2cMapCreate('long', 'int');

for i = 1:nbcs
    bctag = int64(bcpriority(i));
    if m2cMapIsKey(map, bctag)
        map = m2cMapDestroy(map);
        error('build_bcpriority:dupBC', 'duplicated boundary type %d', int32(bctag));
    end
    map = m2cMapSetValue(map, bctag, int32(0));
end
for i = 1:nfacets
    bctag = int64(bctags(i));
    [v, flag] = m2cMapGetValue(map, bctag);
    if flag
        map = m2cMapDestroy(map);
        error('build_bcpriority:unknownBC', 'unknown boundary type %d', int32(bctag));
    end
    map = m2cMapSetValue(map, bctag, v + 1);
end

% build bc list
bc_ptr = m2cNullcopy(zeros(nbcs + 1, 1, 'int32')); bc_ptr(1) = 1;
for i = 1:nbcs
    bctag = int64(bcpriority(i));
    v = m2cMapGetValue(map, bctag);
    if v == 0
        map = m2cMapDestroy(map);
        error('build_bcpriority:emptyBC', 'empty boundary type %d in priority list', int32(bctag));
    end
    bc_ptr(i + 1) = bc_ptr(i) + v;
end
m2cAssert(bc_ptr(end) - 1 == nfacets);
bc_fsets = m2cNullcopy(zeros(bc_ptr(end) - 1, 1, 'int32'));

% reuse map to build bc->bcidx map
for i = 1:nbcs
    bctag = int64(bcpriority(i));
    map = m2cMapSetValue(map, bctag, i);
end

for i = 1:nfacets
    bcidx = m2cMapGetValue(map, int64(bctags(i)));
    bc_fsets(bc_ptr(bcidx)) = fsetids(i);
    bc_ptr(bcidx) = bc_ptr(bcidx) + 1;
end
for i = nbcs:-1:2; bc_ptr(i) = bc_ptr(i - 1); end
bc_ptr(1) = 1;

map = m2cMapDestroy(map);

end

function [mesh, nsetids] = build_nsets(mesh, bc_ptr, bc_list, bcs)
% compute node sets

coder.inline('never');
m2cNowarnBuf('\w+'); % ignore string buffers

nbcs = m2cIgnoreRange(cast(size(bc_ptr, 1) - 1, 'int32'));
nsetids = m2cNullcopy(zeros(nbcs, 1, 'int32'));

uset = m2cSetCreate('int', 'ordered'); % alloc 1

for i = 1:nbcs
    seti = m2cSetCreate('int', 'ordered'); % alloc 2
    for j = bc_ptr(i):bc_ptr(i + 1) - 1
        fsetidx = bc_list(j);
        nfacets = sfemesh_nfacetsets(mesh, fsetidx);
        for k = 1:nfacets
            [eid, facet] = sfemesh_facetset(mesh, fsetidx, k);
            etype = sfemesh_etype(mesh, eid);
            [~, lids] = obtain_facets(etype, facet);
            % local IDs of the element of this facet stored in `lids`
            nn = cast(size(lids, 1), 'int32');
            for ii = 1:nn
                nid = sfemesh_econn(mesh, eid, int32(lids(ii)));
                seti = m2cSetInsert(seti, nid);
            end
        end
    end
    setd = m2cSetDifference(seti, uset); % alloc 3
    [mesh, nsetid] = sfemesh_append_nodeset(mesh);
    nsetids(i) = nsetid;
    mesh.nodesets(nsetid).name = sprintf('BC%d', bcs(i));
    mesh = sfemesh_resize_nodesetdata(mesh, nsetid, m2cSetSize(setd));
    mesh.nodesets(nsetid).nids = m2cSetGetKeys(setd);
    % union
    uset2 = m2cSetUnion(seti, uset); % alloc 4
    uset = m2cSetDestroy(uset); % dealloc 4
    uset = uset2;
    setd = m2cSetDestroy(setd); % deallo 3
    seti = m2cSetDestroy(seti); % dealloc 2
end

uset = m2cSetDestroy(uset); % dealloc 1

end
