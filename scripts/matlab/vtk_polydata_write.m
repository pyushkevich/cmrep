function vtk_polydata_write(file, p)
% Write polydata to a VTK file (ASCII for compatibility)
% Usage:
%   vtk_polydata_write(file, p)

    % open the file
    fid = fopen(file,'w');
    
    % Write the header
    fprintf(fid, '# vtk DataFile Version 3.0\n');
    fprintf(fid, '%s\nASCII\nDATASET POLYDATA\n', p.hdr.name);
    
    % Write the points
    fprintf(fid, 'POINTS %d float\n', size(p.points,1));
    fprintf(fid, '%f %f %f %f %f %f %f %f %f\n', p.points');
    fprintf(fid, '\n');
    
    % Write all the cells
    cell_types = fieldnames(p.cells);
    total_cells = 0;
    for i = 1:length(cell_types)
        pc = p.cells.(cell_types{i});
        n_cell = length(pc);
	    total_cells = total_cells + n_cell;
        n_store = sum(arrayfun(@(x)length(x{1}), pc)) + n_cell;
        fprintf(fid, '%s %d %d\n', cell_types{i}, n_cell, n_store);
        for j = 1:n_cell
            fprintf(fid, '%d', length(pc{j}));
            fprintf(fid, ' %d', pc{j}-1);
            fprintf(fid, '\n');
        end
    end
    
    % Write all the point_data and cell_data attributes
    attr_types = {'point_data', 'cell_data'};
    totals = [size(p.points,1) total_cells];
    for i = 1:length(attr_types)
        
        % Write the entry
        fprintf(fid, '%s %d\n', attr_types{i}, totals(i));
        
        % Get the appropriate parent struct
        if ~isfield(p, attr_types{i}), continue; end
        pa = p.(attr_types{i});
        
        % Write all the field entries
        fld = strcmpi('field',{pa(:).type});
        if any(fld)
            
            % Write the field data header
            fprintf(fid, 'FIELD FieldData %d\n', sum(fld));
            
            % Write each of the fields
            for k = find(fld)                
                fprintf(fid, '%s %d %d double\n', ...
                    vtk_encode(pa(k).name), size(pa(k).data'));
                fprintf(fid, '%f %f %f %f %f %f %f %f %f\n', pa(k).data');
                fprintf(fid, '\n');
            end
            
        end
        
        % Write all the non-field entries
        if ~all(fld)
            
            % Write each of the non-fields
            for k = find(~fld)               
                
                fprintf(fid, '%s', pa(k).type);                
                
                if any(strcmpi(pa(k).type, {'normals','vectors','tensors'}))                    
                    fprintf(fid, ' %s double\n', ...
                        vtk_encode(pa(k).name));
                elseif any(strcmpi(pa(k).type, {'texture_coordinates'})) 
                    fprintf(fid, ' %s %d double\n', pa(k).name, size(pa(k).data,2));
                elseif strcmpi(pa(k).type, 'scalars')
                    fprintf(fid, ' %s float\nlookup_table default\n',pa(k).name); 
                else
                    error('Can not write type %s', pa(k).type);                
                end
                
                fprintf(fid, '%f %f %f %f %f %f %f %f %f\n', pa(k).data');                
                fprintf(fid, '\n');
            end
        end
    
    end
    
    fclose(fid);
    
end

function s = vtk_encode(t)
    s = strrep(t,' ','%20');
end
    
