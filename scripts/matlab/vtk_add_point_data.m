function p1 = vtk_add_point_data(p, name, data, overwrite)
% Add a point attribute array to a vtk mesh
% Usage:
%   p1 = vtk_add_field_data(p, name, data)
% Parameters
%   p         VTK mesh struct (from vtk_polydata_read)
%   name      Name of the new array (string)
%   data      An Nxk matrix of values to add
%   overwrite boolean, whether to overwrite existing array of the same name

arr.name = name;
arr.type = 'field';

if size(data, 1) == size(p.points,1)
    arr.data = data;
elseif size(data, 2) == size(p.points,1)
    arr.data = data';
else
    error('Data size does not match point array size');
end

p1 = p;
if ~isfield(p1, 'point_data')
    p1.point_data(1) = arr;
else
    pos = strmatch(name, {p1.point_data.name}, 'exact');
    if ~isempty(pos)
       if overwrite
          p1.point_data(pos(1)) = arr;
       else
          error('This array already exists');
       end
    else         
       p1.point_data(length(p1.point_data)+1) = arr;
    end
end
