function arr = vtk_get_point_data(p, name)
% Get the point_data array by a certain name from a structure
% returned by vtk_read_polydata
% Usage:
%       arr = vtk_get_point_data(p, name)
index = find(strcmpi(name, {p.point_data(:).name}));
if length(index) ~= 1
    error('Point data not found');
end
arr = p.point_data(index).data;

