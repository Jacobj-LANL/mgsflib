function [edges, lids] = update_edges_table(funcname)
% This script computes the edge mapping with CGNS numbering convention

if nargin < 1 || isempty(funcname); funcname = 'obtain_edges_tet'; end

rt = sfelib_root;

edges = containers.Map('KeyType', 'int32', 'ValueType', 'any');
lids = containers.Map('KeyType', 'int32', 'ValueType', 'any');
names = containers.Map('KeyType', 'int32', 'ValueType', 'char');

[edges, lids, names] = ...
    fetch_data(dir(fullfile(rt, 'sfe/elemtypes')), edges, lids, names);

fid = fopen(fullfile(mgsflib_root, 'demomeshes', [funcname '.m']), 'w');

fprintf(fid, 'function [ret, lids] = %s(etype, edgeid)\n', funcname);
fprintf(fid, '%%%s - Query edge information\n', funcname);
fprintf(fid, '%%\n%%    nedges = %s(etype)\n', funcname);
fprintf(fid, '%%    edgetype = %s(etype, edgeid)\n', funcname);
fprintf(fid, '%%    [edgetype, lids] = %s(etype, edgeid)\n%%\n', funcname);
fprintf(fid, '%% PARAMETERS\n');
fprintf(fid, '%% ----------\n');
fprintf(fid, '%%    etype:          Element type enums\n');
fprintf(fid, '%%    edgeid:        Local edge IDs\n%%\n');
fprintf(fid, '%% RETURNS\n');
fprintf(fid, '%% -------\n');
fprintf(fid, '%% If only etype is provided, then this function returns\n');
fprintf(fid, '%% a single number `nedges`, which is the number\n');
fprintf(fid, '%% of edges of an element.\n%%\n');
fprintf(fid, '%% If edge ID is provided (1:nedges), then this function\n');
fprintf(fid, '%% returns the edge element type of `edgeid`. Optionally,\n');
fprintf(fid, '%% it can also return the local IDs of the requested edge\n%%\n');
fprintf(fid, '%% This file was generated from update_edges_table\n\n');

fprintf(fid, '%%#codegen -args {int32(0)} -nargout 1\n');
fprintf(fid, '%%#codegen -args {int32(0), int8(0)} -nargout 2\n\n');

fprintf(fid, 'persistent EDGES LIDS\n\n');

fprintf(fid, "if coder.target('MATLAB')\n");

% MATLAB code
fprintf(fid, '    if isempty(EDGES); EDGES = gen_edges; end\n');
fprintf(fid, '    if isempty(LIDS); LIDS = gen_lids; end\n');
fprintf(fid, '    if nargin < 2; ret = uint8(size(EDGES(etype), 1)); return; end\n');
fprintf(fid, '    fs = EDGES(etype); ret = fs(edgeid);\n');
fprintf(fid, '    if nargout > 1\n');
fprintf(fid, "        ls = LIDS(etype); lids = int16(ls(edgeid,:)'); lids = lids(lids~=0);\n");
fprintf(fid, '    end\n');
fprintf(fid, '    return;\n');
fprintf(fid, 'end  %% MATLAB code ends\n\n');

% C++ or MEX
fprintf(fid, "coder.cinclude('<vector>');\n");
fprintf(fid, "coder.cinclude('<algorithm>');\n");
fprintf(fid, "coder.inline('never');\n\n");

% edges
fprintf(fid, "coder.ceval('const static std::vector<std::vector<uint8_T>> EDGES{//');\n");

edge_keys = keys(edges);
start = min([edge_keys{:}]);
for i = start:max([edge_keys{:}])
    if ~edges.isKey(i)
        fprintf(fid, "coder.ceval('{},// %d');\n", i);
    else
        fprintf(fid, "coder.ceval('{");
        fs = edges(i);
        for j = 1:numel(fs) - 1
            fprintf(fid, '%d,', fs(j));
        end
        fprintf(fid, "%d},// %s');\n", fs(end), names(i));
    end
end

fprintf(fid, "coder.ceval('};//');\n"); % initializer for edges

% local IDs
fprintf(fid, 'if nargin > 1\n');
fprintf(fid, "    coder.ceval('const static std::vector<std::vector<std::vector<int16_T>>> LIDS{//');\n");

lid_keys = keys(lids);
for i = start:max([lid_keys{:}])
    if ~lids.isKey(i)
        fprintf(fid, "    coder.ceval('{{}},// %d');\n", i);
    else
        fprintf(fid, "    coder.ceval('{");
        ls = lids(i);
        [m, n] = size(ls);
        for j = 1:m
            fprintf(fid, '{');
            for k = 1:n - 1
                fprintf(fid, '%d,', ls(j, k));
            end
            if j < m
                fprintf(fid, '%d},', ls(j, end));
            else
                fprintf(fid, '%d}', ls(j, end));
            end
        end
        fprintf(fid, "},// %s');\n", names(i));
    end
end

fprintf(fid, "    coder.ceval('};//');\n"); % initializer for lids
fprintf(fid, 'end\n');

fprintf(fid, 'n = int32(0);\n');
fprintf(fid, 'ret = uint8(0);\n');

fprintf(fid, 'if nargin < 2\n');
fprintf(fid, '%% get the number of edges\n');
fprintf(fid, "    ret = coder.ceval('[&](uint8_T et){return EDGES[et-%d].size();}', etype); return;\n", start);
fprintf(fid, 'end\n');

fprintf(fid, "ret = coder.ceval('[&](int et, uint8_T fid){return EDGES[et-%d][fid];}', etype, edgeid-1);\n", start);

fprintf(fid, 'if nargout > 1\n');
fprintf(fid, "    n = coder.ceval('[&](int et, uint8_T fid){int n = LIDS[et-%d][fid].size(); while (n && LIDS[et-%d][fid][n-1] == 0) --n; return n;}', etype, edgeid-1);\n", start, start);
fprintf(fid, "    lids = coder.nullcopy(zeros(n, 1, 'int16'));\n");
fprintf(fid, "    coder.varsize('lids', [50 1], [1 0]);\n");
fprintf(fid, "    coder.ceval('[&](int et, uint8_T fid, int n, std::int16_t *v){std::copy_n(LIDS[et-%d][fid].cbegin(), n, v);}', etype, edgeid-1, n, coder.wref(lids));\n", start);
fprintf(fid, 'end\n');

fprintf(fid, '\nend\n');

% local function for initializing edges
fprintf(fid, '\nfunction edges = gen_edges\n');
fprintf(fid, "if coder.target('MATLAB')\n");
fprintf(fid, '    ks = {');
for i = 1:numel(edge_keys) - 1
    fprintf(fid, 'uint8(%d),', edge_keys{i});
end
fprintf(fid, 'uint8(%d)};\n', edge_keys{end});
fprintf(fid, '    vs = {');
vs = values(edges);
for i = 1:numel(vs) - 1
    fprintf(fid, 'uint8(['); fprintf(fid, '%d;', vs{i}); fprintf(fid, ']),...\n        ');
end
fprintf(fid, 'uint8(['); fprintf(fid, '%d;', vs{end}); fprintf(fid, '])};\n');
fprintf(fid, '    edges = containers.Map(ks, vs);\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end\n');

% local function for initializing local IDs
fprintf(fid, '\nfunction lids = gen_lids\n');
fprintf(fid, "if coder.target('MATLAB')\n");
fprintf(fid, '    ks = {');
for i = 1:numel(lid_keys) - 1
    fprintf(fid, 'uint8(%d),', lid_keys{i});
end
fprintf(fid, 'uint8(%d)};\n', lid_keys{end});
fprintf(fid, '    vs = {');
vs = values(lids);
for i = 1:numel(vs) - 1
    fprintf(fid, 'uint8([');
    ls = vs{i};
    [m, n] = size(ls);
    for j = 1:m
        if n > 1
            fprintf(fid, '%d,', ls(j, 1:end - 1));
        end
        fprintf(fid, '%d;', ls(j, end));
    end
    fprintf(fid, ']),...\n        ');
end
ls = vs{end};
[m, n] = size(ls);
fprintf(fid, 'uint8([');
for j = 1:m
    if n > 1
        fprintf(fid, '%d,', ls(j, 1:end - 1));
    end
    fprintf(fid, '%d;', ls(j, end));
end
fprintf(fid, '])};\n');
fprintf(fid, '    lids = containers.Map(ks, vs);\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end\n');

fclose(fid);

clear(funcname);

end

function [data, lids, names] = fetch_data(fs, data, lids, names)

for i = 1:numel(fs)
    if startsWith(fs(i).name, 'SFE_') && ~startsWith(fs(i).name, 'SFE_SHAPE')
        pat = split(fs(i).name, '.');
        fname = pat{end - 1};
        etype = feval(fname);
        if isKey(data, etype); continue; end
        if etype == 1; continue; end % skip for nodes
        dim = obtain_elemdim(etype);
        shape = obtain_elemshape(etype);
        bins = dec2bin(etype, 8);
        deg = obtain_elemdegree(etype);
        names(etype) = fname;
        bintmp = bins;
        if shape == 4
            % tetras
            bintmp(1:3) = '001';
            facetbar = to_equi_if_necessary(bintmp);
            data(etype) = [facetbar; facetbar; facetbar; facetbar; facetbar; facetbar];
            switch deg
                case 1
                    lids(etype) = uint8([1,2; 2,3; 3,1; 1,4; 2,4; 3,4]);
                case 2
                    lids(etype) = uint8([1,2,5; 2,3,6; 3,1,7; 1,4,8; 2,4,9; 3,4,10]);
                case 3
                    lids(etype) = uint8([1,2,5,6; 2,3,7,8; 3,1,9,10; 1,4,11,12; 2,4,13,14; 3,4,15,16]);
                case 4
                    lids(etype) = uint8([1,2,5,6,7; 2,3,8,9,10; 3,1,11,12,13; 1,4,14,15,16; 2,4,17,18,19; 3,4,20,21,22]);
                case 5
                    lids(etype) = uint8([1,2,5,6,7,8; 2,3,9,10,11,12; 3,1,13,14,15,16; 1,4,17,18,19,20; 2,4,21,22,23,24; 3,4,25,26,27,28]);
                case 6
                    lids(etype) = uint8([1,2,5,6,7,8,9; 2,3,10,11,12,13,14; 3,1,15,16,17,18,19; 1,4,20,21,22,23,24; 2,4,25,26,27,28,29; 3,4,30,31,32,33,34]);
                otherwise
                    assert(false, 'not supported');
            end
        end
    end
end

end

function etype = to_equi_if_necessary(bins)
etype = uint8(bin2dec(bins));
if obtain_elemnnodes(etype) == 0
    % use equidistant
    bins(end - 1:end) = '00';
    etype = uint8(bin2dec(bins));
end
assert(obtain_elemnnodes(etype));
end
