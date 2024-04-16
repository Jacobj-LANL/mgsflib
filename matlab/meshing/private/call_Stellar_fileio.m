function [tets,coords] = call_Stellar_fileio(tets, coords,bdy, remake)

if nargin < 4
    remake = false;
end


currdir = pwd;
dir = [mgsflib_root, '/meshing/private/Stellar'];
cd(dir);
if remake
    system('rm -f Stellar');
    system('rm -f Starbase*');
end
if ~isfile([dir,'/Stellar'])
    system('make')
end

% write .ele file
fid = fopen('tets.ele','w');
nelems = size(tets,1);
fprintf(fid, '%d 4 0\n',nelems);
fprintf(fid, '\t%d %d %d %d %d\n',[(1:nelems)', tets(:,[2,3,4,1])]');
fprintf(fid, '#Generated by Call_Stellar_filio.m');
fclose(fid);

%write .node file
fid = fopen('tets.node','w');
nv = size(coords,1);
fprintf(fid, '%d 3 0 1\n',nv);
fprintf(fid, '\t%d %0.16f %0.16f %0.16f %d\n',[(1:nv)', coords, bdy]');
fprintf(fid, '#Generated by Call_Stellar_filio.m');
fclose(fid);

system('./Stellar -s EXAMPLE_CONFIG tets')

% read files
fid = fopen('tets.1.ele', 'r');
A = textscan(fid, '%d %d %d %d %d %f', 'HeaderLines', 1, 'Delimiter', ' ', 'CommentStyle', 'Shell', 'MultipleDelimsAsOne', 1);
fclose(fid);
eids = int32(A{1});
elems = zeros(max(eids), 4, 'int32');
elems(eids, :) = int32([A{2} A{3} A{4} A{5}]);
tets = elems(:,[2,3,4,1]);

fid = fopen('tets.1.node', 'r');
[A] = textscan(fid, '%d %f %f %f %d', 'HeaderLines', 1, 'Delimiter', ' ', 'CommentStyle', 'Shell', 'MultipleDelimsAsOne', 1);
fclose(fid);
nids = double(A{1});
xs = zeros(max(nids), 3);
xs(nids, :) = [A{2} A{3} A{4}];
coords = xs;
cd(currdir)
end
