function codegen_mgsflib(varargin)
% codegen_mgsflib - Generate C++ libraries for submodules
%
%     codegen_mgsflib <codegen_lib_args> [module]
%
% Examples:
%     codegen_mgsflib -force hiermesh
%     codegen_mgsflib -force sfe
%     codegen_mgsflib -force sfe/shapefuncs
%     codegen_mgsflib -force
%     codegen_mgsflib -j
%
% Notes
% ----
% Only files in the select subfolders that contain `%#codegen -args` or
% `%#codegen -mex -args` and used by other modules will be compiled.
% For fast compilation, each submodule is compiled into a single MEX
% file named `private/<module>_gateway.mex*` using pre-generated C++
% code. To regenerated C++ code, pass the `-force` or `-j` option.
% The latter will run multiple codegen instances in parallel.
%
% See also build_mgsflib

curpath = pwd;
cleanup = onCleanup(@()cd(curpath));
cd(mgsflib_root);

if nargin && varargin{end}(1) ~= '-'
    submodules = varargin(end);
    hasmodule = 1;
    assert(strcmp(submodules, 'meshing') || strcmp(submodules, 'sfe_assembly') || strcmp(submodules, 'aes_assembly') ...
        || strcmp(submodules, 'stencils'));
else
    submodules = {'meshing','sfe_assembly','aes_assembly','stencils'};
    hasmodule = 0;
end

incs = mgsflib_includes;
nthreads = int32(any(strcmp(varargin, '-j'))) * min(ompGetMaxThreads / 2, 8);
parfor (i = 1:length(submodules), nthreads)
    folder = [submodules{i} '/'];
    files = grep_files([folder '*.m'], '\n%#codegen\s+(-mex\s+)?-args');

    [~, basename] = fileparts(submodules{i});

    codegen_lib('-mex', '-O3', '-replace-types', '-no-report', incs{:}, ...
        varargin{1:end - hasmodule}, [folder 'private/' basename '_gateway.m'], files{:}); %#ok<PFBNS>
end

end