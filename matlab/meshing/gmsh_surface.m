function [tris, xs, normals, inside] = gmsh_surface(modelname, n, min_sz)
persistent compiled;
dir = [ mgsflib_root '/meshing'];
if ~exist([dir '/private/gmsh_mex.' mexext], 'file') || (nargin==0 || isempty(compiled)) || ...
        isnewer([dir '/private/gmsh_install/gmsh_mex.cpp'],which('gmsh_mex'))
    disp('Building gmsh_surface...');
    command = ['mex -g -I',dir,'/private/gmsh_install/include ', '-outdir ',dir,'/private ', dir, '/private/gmsh_install/gmsh_mex.cpp ','-L. ',dir,'/private/gmsh_install/lib/libgmsh.a'];
    eval(command);
    compiled = 1;
    disp('gmsh_surface was built successfully.');
end

if nargin == 0
    return;
end

if nargout == 3
    [tris,xs,normals] = gmsh_mex(n,modelname,min_sz);
else
    [tris,xs,normals,inside] = gmsh_mex(n,modelname,min_sz);
end

end
