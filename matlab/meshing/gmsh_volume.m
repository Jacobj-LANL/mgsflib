function [tets, xs] = gmsh_volume(modelname, n, min_sz)

dir = [ mgsflib_root '/meshing'];
if ~exist([dir '/private/gmsh_tetmex.' mexext], 'file') || (nargin==0) || ...
        isnewer([dir '/private/gmsh_install/gmsh_tetmex.cpp'],which('gmsh_tetmex'))
    disp('Building gmsh_volume...');
    command = ['mex -g -I',dir,'/private/gmsh_install/include ', '-outdir ',dir,'/private ', dir, '/private/gmsh_install/gmsh_tetmex.cpp ','-L. ',dir,'/private/gmsh_install/lib/libgmsh.a'];
    eval(command);
    compiled = 1;
    disp('gmsh_volume was built successfully.');
end

if nargin == 0
    return;
end

[tets,xs] = gmsh_tetmex(n,modelname,min_sz);
end