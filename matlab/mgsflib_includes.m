function incs = mgsflib_includes
% Returns a cell array of include paths
incs = {['-I"' fullfile(fileparts(sfelib_root), 'cpp', 'src') '"'], ...
    ['-I"' fullfile(fileparts(wlslib_root), 'cpp', 'src') '"']};