function p = vtk_polydata_read(file, varargin)
% Read VTK polydata into a struct
% Usage:
%   p = vtk_polydata_read(file, pars)
% Parameters:
%   pars.encoding           One of 'ieee-le', 'ieee-be' (default)
% Return Value:
%   p.hdr                   Header information
%   p.points                N x 3 array of point coordinates
%   p.cells.(name)          cell array of polygons, lines, etc
%   p.point_data            Point data arrays, following format:
%   p.point_data(i).name    Array name
%   p.point_data(i).type    VTK type (normals,  vectors, field)
%   p.point_data(i).data    N x K array of values

    % Initialize the parameters
    if(nargin > 1) pars = varargin{1}; else pars = struct(); end;
    if(~isfield(pars,'encoding')) pars.encoding='ieee-be'; end;        

    % Open file
    fid = fopen(file, 'r');
    
    % Read the header (2 lines)
    p.hdr.name = vtkreadstr(fid, '%[^\n]\n');
    p.hdr.type = vtkreadstr(fid, '%s\n');
      
    % Read the dataset type
    p.hdr.dst = vtkreadstr(fid, 'DATASET %s\n');
    if ~strcmpi(p.hdr.dst, 'POLYDATA')
        error('Data types other than POLYDATA unsupported');
    end
    
    % Read the rest of the stuff in the file
    mode = 'intro';
    nattr = 0;
    
    while ~feof(fid)
        
        % Get the keyword (first word of the line)
        key = vtkreadstr(fid, '%s');        
        
        if (length(key) == 0) continue; end
        
        % fprintf('Reading %s\n', vtk_decode(key));                    
            
        if strcmp(mode,'intro') && strcmpi(key, 'points')
            
            % Read the number of points and the type
            dat = vtkread(fid, '%d %s\n', 1);
            n = double(dat{1});
            data_type = char(dat{2});

            % Read point values
            X = vtkreaddata(fid, p, 3 * n, data_type, pars);
            p.points = reshape(X, 3, [])';

        elseif strcmp(mode,'intro') && any(strcmpi(key, ...
                {'polygons','vertices','lines','triangle_strips'}))
                
            % Read the number of cells and storage size
            dat = vtkread(fid, '%d %d\n', 1);
            n = double(dat{1});
            storage = double(dat{2});
            data_type = 'uint';            
            
            % Read the cell data
            T = vtkreaddata(fid, p, storage, data_type, pars);
            
            % Connect the cell data together
            i = 1; j = 1;
            while i <= storage
                cd{j} = T(i+1: i+T(i)) + 1;
                j = j + 1;
                i = i + T(i) + 1;
            end
               
            % Place in the appropriate array
            p.cells.(lower(key)) = cd;
            
        elseif any(strcmpi(key, {'point_data', 'cell_data'}))
            
            % Set the mode variable
            mode = lower(key);
            
            % Initialize the data arrays
            % p.(mode) = [];
            
            % Reset the field counter
            n_fields = 0;
            
            % Read the number of attributes
            nattr = vtkreadnum(fid, '%d\n');
                        
        elseif any(strcmpi(mode, {'point_data', 'cell_data'}))
            
            if any(strcmpi(key, {'normals','scalars','color_scalars',...
                'vectors', 'texture_coordinates', 'tensors'}))
            
                % Create an array
                arr.name = vtk_decode(vtkreadstr(fid, '%s'));
                arr.type = lower(key);

                % We are reading some sort of attributes
                if any(strcmpi(key, {'normals','vectors'}))

                    % Read the data
                    data_type = vtkreadstr(fid, '%s\n');
                    X = vtkreaddata(fid, p, nattr * 3, data_type, pars);
                    arr.data = reshape(X, 3, [])';

                elseif strcmpi(key, 'tensors')

                    % Read the data
                    data_type = vtkreadstr(fid, '%s\n');
                    X = vtkreaddata(fid, p, nattr * 9, data_type, pars);
                    arr.data = reshape(X, 9, [])';

                elseif strcmpi(key, 'texture_coordinates')

                    % Read the data
                    ncomp = vtkreadnum(fid, '%d');
                    data_type = vtkreadstr(fid, '%s\n');
                    X = vtkreaddata(fid, p, nattr * ncomp, data_type, pars);
                    arr.data = reshape(X, ncomp, [])';

                elseif strcmpi(key, 'scalars')
                    
                    % Read the scarar data
                    data_type = vtkreadstr(fid, '%s\n');
                    junk = vtkreadstr(fid, '%s %s\n');
                    X = vtkreaddata(fid, p, nattr, data_type, pars);
                    arr.data = reshape(X, 1, [])';
                    
                else
                    
                    error('Dataformat %s is not yet supported', key);

                end

                % Append the array
                if ~isfield(p, mode)
                    p.(mode)(1) = arr;
                else
                    p.(mode)(1+length(p.(mode))) = arr;
                end
            
            elseif strcmpi(key, 'field')
                
                % Enter field reading mode
                vtkreadstr(fid, '%s');
                n_fields = vtkreadnum(fid, '%d\n');
                
            elseif n_fields > 0
                
                % Read the field information
                arr.name = vtk_decode(key);
                arr.type = 'field';
                ncomp = vtkreadnum(fid, '%d');
                ntuples = vtkreadnum(fid, '%d');
                data_type = vtkreadstr(fid, '%s\n');
                
                % Read the tuple data
                X = vtkreaddata(fid, p, ncomp * ntuples, data_type, pars);
                arr.data = reshape(X, ncomp, [])';
                                
                % Append the array
                if ~isfield(p, mode)
                    p.(mode)(1) = arr;
                else
                    p.(mode)(1+length(p.(mode))) = arr;
                end
                
                % Decrement the field counter
                n_fields = n_fields - 1;
                
            else
                
                error('Unknown entry %s', key);
                
            end
            
        else
            
            error('Unknown entry %s', key);

        end
    end


    fclose(fid);
end

function s = vtkread(fid, str, n)
    s = textscan(fid, str, n, 'ReturnOnError', 0, 'CommentStyle', '#');
end

function s = vtkreadstr(fid, pat)
    str = textscan(fid, pat, 1, 'ReturnOnError', 0, 'CommentStyle', '#');
    s = char(str{1});
end

function n = vtkreadnum(fid, pat)
    val = textscan(fid, pat, 1, 'ReturnOnError', 0, 'CommentStyle', '#');
    n = double(val{1});
end

function X = vtkreaddata(fid, p, comp, data_type, pars)
    if strcmpi(p.hdr.type, 'ascii')
        car = textscan(fid, '%f', comp, ...
            'ReturnOnError', 0, 'CommentStyle', '#');
        X = car{1};
    else
        if strcmpi(data_type, 'float')
            conv = 'float32=>double';
        elseif strcmpi(data_type, 'double')
            conv = 'float64=>double';
        elseif strcmpi(data_type, 'uint')
            conv = 'uint32=>uint32';
        elseif strcmpi(data_type, 'vtkIdType')
            conv = 'uint32=>uint32';
        else
            error('data_type %s is unsupported', data_type);
        end
            
        X = fread(fid, comp, conv, pars.encoding);
    end
end
    
function S = vtk_decode(t)
    S = strrep(t,'%20',' ');
end
