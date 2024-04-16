function [ret, lids] = obtain_edges_tet(etype, edgeid)
%obtain_edges_tet - Query edge information
%
%    nedges = obtain_edges_tet(etype)
%    edgetype = obtain_edges_tet(etype, edgeid)
%    [edgetype, lids] = obtain_edges_tet(etype, edgeid)
%
% PARAMETERS
% ----------
%    etype:          Element type enums
%    edgeid:        Local edge IDs
%
% RETURNS
% -------
% If only etype is provided, then this function returns
% a single number `nedges`, which is the number
% of edges of an element.
%
% If edge ID is provided (1:nedges), then this function
% returns the edge element type of `edgeid`. Optionally,
% it can also return the local IDs of the requested edge
%
% This file was generated from update_edges_table

%#codegen -args {int32(0)} -nargout 1
%#codegen -args {int32(0), int8(0)} -nargout 2

persistent EDGES LIDS

if coder.target('MATLAB')
    if isempty(EDGES); EDGES = gen_edges; end
    if isempty(LIDS); LIDS = gen_lids; end
    if nargin < 2; ret = uint8(size(EDGES(etype), 1)); return; end
    fs = EDGES(etype); ret = fs(edgeid);
    if nargout > 1
        ls = LIDS(etype); lids = int16(ls(edgeid,:)'); lids = lids(lids~=0);
    end
    return;
end  % MATLAB code ends

coder.cinclude('<vector>');
coder.cinclude('<algorithm>');
coder.inline('never');

coder.ceval('const static std::vector<std::vector<uint8_T>> EDGES{//');
coder.ceval('{36,36,36,36,36,36},// SFE_TET_4');
coder.ceval('{},// 133');
coder.ceval('{},// 134');
coder.ceval('{},// 135');
coder.ceval('{40,40,40,40,40,40},// SFE_TET_10');
coder.ceval('{},// 137');
coder.ceval('{},// 138');
coder.ceval('{},// 139');
coder.ceval('{44,44,44,44,44,44},// SFE_TET_20');
coder.ceval('{45,45,45,45,45,45},// SFE_TET_FEK_20');
coder.ceval('{},// 142');
coder.ceval('{},// 143');
coder.ceval('{48,48,48,48,48,48},// SFE_TET_35');
coder.ceval('{49,49,49,49,49,49},// SFE_TET_GL_35');
coder.ceval('{48,48,48,48,48,48},// SFE_TET_FEK_35');
coder.ceval('{},// 147');
coder.ceval('{52,52,52,52,52,52},// SFE_TET_56');
coder.ceval('{53,53,53,53,53,53},// SFE_TET_GL_56');
coder.ceval('{52,52,52,52,52,52},// SFE_TET_FEK_56');
coder.ceval('{},// 151');
coder.ceval('{56,56,56,56,56,56},// SFE_TET_84');
coder.ceval('{57,57,57,57,57,57},// SFE_TET_GL_84');
coder.ceval('{56,56,56,56,56,56},// SFE_TET_FEK_84');
coder.ceval('};//');
if nargin > 1
    coder.ceval('const static std::vector<std::vector<std::vector<int16_T>>> LIDS{//');
    coder.ceval('{{1,2},{2,3},{3,1},{1,4},{2,4},{3,4}},// SFE_TET_4');
    coder.ceval('{{}},// 133');
    coder.ceval('{{}},// 134');
    coder.ceval('{{}},// 135');
    coder.ceval('{{1,2,5},{2,3,6},{3,1,7},{1,4,8},{2,4,9},{3,4,10}},// SFE_TET_10');
    coder.ceval('{{}},// 137');
    coder.ceval('{{}},// 138');
    coder.ceval('{{}},// 139');
    coder.ceval('{{1,2,5,6},{2,3,7,8},{3,1,9,10},{1,4,11,12},{2,4,13,14},{3,4,15,16}},// SFE_TET_20');
    coder.ceval('{{1,2,5,6},{2,3,7,8},{3,1,9,10},{1,4,11,12},{2,4,13,14},{3,4,15,16}},// SFE_TET_FEK_20');
    coder.ceval('{{}},// 142');
    coder.ceval('{{}},// 143');
    coder.ceval('{{1,2,5,6,7},{2,3,8,9,10},{3,1,11,12,13},{1,4,14,15,16},{2,4,17,18,19},{3,4,20,21,22}},// SFE_TET_35');
    coder.ceval('{{1,2,5,6,7},{2,3,8,9,10},{3,1,11,12,13},{1,4,14,15,16},{2,4,17,18,19},{3,4,20,21,22}},// SFE_TET_GL_35');
    coder.ceval('{{1,2,5,6,7},{2,3,8,9,10},{3,1,11,12,13},{1,4,14,15,16},{2,4,17,18,19},{3,4,20,21,22}},// SFE_TET_FEK_35');
    coder.ceval('{{}},// 147');
    coder.ceval('{{1,2,5,6,7,8},{2,3,9,10,11,12},{3,1,13,14,15,16},{1,4,17,18,19,20},{2,4,21,22,23,24},{3,4,25,26,27,28}},// SFE_TET_56');
    coder.ceval('{{1,2,5,6,7,8},{2,3,9,10,11,12},{3,1,13,14,15,16},{1,4,17,18,19,20},{2,4,21,22,23,24},{3,4,25,26,27,28}},// SFE_TET_GL_56');
    coder.ceval('{{1,2,5,6,7,8},{2,3,9,10,11,12},{3,1,13,14,15,16},{1,4,17,18,19,20},{2,4,21,22,23,24},{3,4,25,26,27,28}},// SFE_TET_FEK_56');
    coder.ceval('{{}},// 151');
    coder.ceval('{{1,2,5,6,7,8,9},{2,3,10,11,12,13,14},{3,1,15,16,17,18,19},{1,4,20,21,22,23,24},{2,4,25,26,27,28,29},{3,4,30,31,32,33,34}},// SFE_TET_84');
    coder.ceval('{{1,2,5,6,7,8,9},{2,3,10,11,12,13,14},{3,1,15,16,17,18,19},{1,4,20,21,22,23,24},{2,4,25,26,27,28,29},{3,4,30,31,32,33,34}},// SFE_TET_GL_84');
    coder.ceval('{{1,2,5,6,7,8,9},{2,3,10,11,12,13,14},{3,1,15,16,17,18,19},{1,4,20,21,22,23,24},{2,4,25,26,27,28,29},{3,4,30,31,32,33,34}},// SFE_TET_FEK_84');
    coder.ceval('};//');
end
n = int32(0);
ret = uint8(0);
if nargin < 2
% get the number of edges
    ret = coder.ceval('[&](uint8_T et){return EDGES[et-132].size();}', etype); return;
end
ret = coder.ceval('[&](int et, uint8_T fid){return EDGES[et-132][fid];}', etype, edgeid-1);
if nargout > 1
    n = coder.ceval('[&](int et, uint8_T fid){int n = LIDS[et-132][fid].size(); while (n && LIDS[et-132][fid][n-1] == 0) --n; return n;}', etype, edgeid-1);
    lids = coder.nullcopy(zeros(n, 1, 'int16'));
    coder.varsize('lids', [50 1], [1 0]);
    coder.ceval('[&](int et, uint8_T fid, int n, std::int16_t *v){std::copy_n(LIDS[et-132][fid].cbegin(), n, v);}', etype, edgeid-1, n, coder.wref(lids));
end

end

function edges = gen_edges
if coder.target('MATLAB')
    ks = {uint8(132),uint8(136),uint8(140),uint8(141),uint8(144),uint8(145),uint8(146),uint8(148),uint8(149),uint8(150),uint8(152),uint8(153),uint8(154)};
    vs = {uint8([36;36;36;36;36;36;]),...
        uint8([40;40;40;40;40;40;]),...
        uint8([44;44;44;44;44;44;]),...
        uint8([45;45;45;45;45;45;]),...
        uint8([48;48;48;48;48;48;]),...
        uint8([49;49;49;49;49;49;]),...
        uint8([48;48;48;48;48;48;]),...
        uint8([52;52;52;52;52;52;]),...
        uint8([53;53;53;53;53;53;]),...
        uint8([52;52;52;52;52;52;]),...
        uint8([56;56;56;56;56;56;]),...
        uint8([57;57;57;57;57;57;]),...
        uint8([56;56;56;56;56;56;])};
    edges = containers.Map(ks, vs);
end
end

function lids = gen_lids
if coder.target('MATLAB')
    ks = {uint8(132),uint8(136),uint8(140),uint8(141),uint8(144),uint8(145),uint8(146),uint8(148),uint8(149),uint8(150),uint8(152),uint8(153),uint8(154)};
    vs = {uint8([1,2;2,3;3,1;1,4;2,4;3,4;]),...
        uint8([1,2,5;2,3,6;3,1,7;1,4,8;2,4,9;3,4,10;]),...
        uint8([1,2,5,6;2,3,7,8;3,1,9,10;1,4,11,12;2,4,13,14;3,4,15,16;]),...
        uint8([1,2,5,6;2,3,7,8;3,1,9,10;1,4,11,12;2,4,13,14;3,4,15,16;]),...
        uint8([1,2,5,6,7;2,3,8,9,10;3,1,11,12,13;1,4,14,15,16;2,4,17,18,19;3,4,20,21,22;]),...
        uint8([1,2,5,6,7;2,3,8,9,10;3,1,11,12,13;1,4,14,15,16;2,4,17,18,19;3,4,20,21,22;]),...
        uint8([1,2,5,6,7;2,3,8,9,10;3,1,11,12,13;1,4,14,15,16;2,4,17,18,19;3,4,20,21,22;]),...
        uint8([1,2,5,6,7,8;2,3,9,10,11,12;3,1,13,14,15,16;1,4,17,18,19,20;2,4,21,22,23,24;3,4,25,26,27,28;]),...
        uint8([1,2,5,6,7,8;2,3,9,10,11,12;3,1,13,14,15,16;1,4,17,18,19,20;2,4,21,22,23,24;3,4,25,26,27,28;]),...
        uint8([1,2,5,6,7,8;2,3,9,10,11,12;3,1,13,14,15,16;1,4,17,18,19,20;2,4,21,22,23,24;3,4,25,26,27,28;]),...
        uint8([1,2,5,6,7,8,9;2,3,10,11,12,13,14;3,1,15,16,17,18,19;1,4,20,21,22,23,24;2,4,25,26,27,28,29;3,4,30,31,32,33,34;]),...
        uint8([1,2,5,6,7,8,9;2,3,10,11,12,13,14;3,1,15,16,17,18,19;1,4,20,21,22,23,24;2,4,25,26,27,28,29;3,4,30,31,32,33,34;]),...
        uint8([1,2,5,6,7,8,9;2,3,10,11,12,13,14;3,1,15,16,17,18,19;1,4,20,21,22,23,24;2,4,25,26,27,28,29;3,4,30,31,32,33,34;])};
    lids = containers.Map(ks, vs);
end
end
