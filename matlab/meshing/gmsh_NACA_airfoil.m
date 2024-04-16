function gmsh_NACA_airfoil(digits, scale_pts, extrude_dim, np)
% this function creates a geo model for the NACA airfoil of specified
% digits and saves file

if nargin < 4
    np = 200;
end
if nargin < 3
    extrude_dim = 3; % which dimension the model is extruded in, 1 for x, 2 for y, and 3 for z
end
if nargin < 2
    scale_pts = [0.45,0.55; 0.4, 0.6]; % first range is the domain of airfoil x
    % second range is the extrusion
end
if nargin < 1
    digits = '2412'; % airfoil code
end

foil_pts = NACA_airfoil(digits, [scale_pts(1,1:2),0.5],np);
nv = size(foil_pts,1);

coords = zeros(nv,3);
dims = setdiff([1,2,3],extrude_dim);
coords(:,dims(1)) = foil_pts(:,1);
coords(:,dims(2)) = foil_pts(:,2);
coords(:,extrude_dim) = scale_pts(2,1);

filename = [mgsflib_root,'/meshing/private/models/NACA_airfoil_',digits,'.geo'];
fid = fopen(filename,'w');

fprintf(fid,'cl__1 = %g;\n',(scale_pts(1,2)-scale_pts(1,1))/50);
for ii = 1:nv
    fprintf(fid,'Point(%d) = {%g, %g, %g, cl__1};\n',ii,coords(ii,1),coords(ii,2),coords(ii,3));
end

fprintf(fid,'Spline(1) = {');
for ii = 1:np-1
    fprintf(fid,'%d,',ii);
end
fprintf(fid,'%d};\n',np);

fprintf(fid,'Spline(2) = {');
for ii = np:nv
    fprintf(fid,'%d,',ii);
end
fprintf(fid,'%d};\n',1);

fprintf(fid, 'Curve Loop(1) = {1,2};\n');
fprintf(fid, 'Plane Surface(1) = {1};\n');

extrude = zeros(1,3);
extrude(extrude_dim) = scale_pts(2,2)-scale_pts(2,1);
fprintf(fid, 'Extrude {%g, %g, %g} {\n\tSurface{1};\n}',extrude(1),extrude(2),extrude(3));
fclose(fid);
end