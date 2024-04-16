function build_mgsflib(varargin)
% build_sfelib - Generate C++ libraries for submodules
%
%     build_sfelib [-force|-j] [module]
%
% Examples:
%     build_sfelib
%     build_sfelib meshdb
%     build_sfelib -force
%     build_sfelib -j
%
% Notes
% ----
% Only files in the select subfolders that contain `%#codegen -args` or
% `%#codegen -mex -args` and used by other modules will be compiled.
% For fast compilation, each submodule is compiled into a single MEX
% file named `private/<module>_gateway.mex*` using pre-generated C++
% code. To regenerated C++ code, run `codegen_sfelib -force [module]`.
%
% See also codegen_sfelib

curpath = pwd;
cleanup = onCleanup(@()cd(curpath));
cd(mgsflib_root);

if nargin && varargin{end}(1) ~= '-'
    submodules = varargin(end);
    assert(strcmp(submodules, 'meshing') || strcmp(submodules, 'sfe_assembly') || strcmp(submodules, 'aes_assembly'));
else
    submodules = {'meshing','sfe_assembly','aes_assembly'};
end

nthreads = int32(any(strcmp(varargin, '-j'))) * min(ompGetMaxThreads / 2, 8);
parfor (i = 1:length(submodules), nthreads)
    folder = [submodules{i} '/'];
    modname = submodules{i};
    [~, basename] = fileparts(modname);
    if any(strcmp(varargin, '-force')) || any(strcmp(varargin, '-j')) || ...
            isnewer(fullfile([folder, 'private'], 'codegen', 'lib', [basename '_gateway'], ['mex_' basename '_gateway.m']), ...
            fullfile([folder, 'private'], [basename '_gateway.', mexext]))
        runscript(fullfile([folder, 'private'], 'codegen', 'lib', [basename '_gateway'], ['mex_' basename '_gateway.m']));
    end
end

end

function runscript(file)
run(file);
end
