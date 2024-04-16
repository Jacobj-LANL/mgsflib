function rt = mgsflib_root
% SFELIB root folder

persistent root__;

if isempty(root__)
    root__ = fileparts(which('mgsflib_root'));
end

rt = root__;

end