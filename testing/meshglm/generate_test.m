%% Load the surface data
m = cell(2,1);
m{1} = vtk_polydata_read('pmatlas_cmrep.bnd.vtk');
m{2} = vtk_mesh_read('pmatlas_cmrep_tet.vtk');
n = size(m{1}.points,1);
T = cell2mat(m{1}.cells.polygons);

%% Generate random collection of 40 subjects
stream = RandStream.getDefaultStream();
stream.reset(2222);

ns = 40;
X = zeros(ns,3);
X(1:20,1) = 1;          % 20 subjects in group 1
X(21:40,2) = 1;         % 20 subjects in group 2
X(:,3) = rand(40,1);    % Third variable random between 0 and 1
con = [-1 1 0.05];

%% Select a set of control points near the mesh. 
% These are randomly chosen points in the dataset, with randomly 
% chosen radius values. They are used to construct a random field 
% over the rest of the image
ncp = 20;
xcp = m{1}.points(randi(n,[ncp ,1]),:) + randn(ncp, 3);
rcp = rand(ncp, 1) * 5;
rcp(1) = 4;
rcp(2) = 2;

% Generate random data at the control points. We let two control points
% be dependent on the subject data, while the others are random
sigma = 2;
fcp = sigma * randn(ncp, ns);
fcp(1,:) = fcp(1,:) + con * X';
fcp(2,:) = fcp(2,:) + con * X';

tcp = zeros(ncp, 1);
tcp(1,:) = 1;
tcp(2,:) = 1;

%% Interpolate f and the true contrast at each point
m1 = cell(2,1);
for q = 1:2
    
    % Interpolate f and the truth at each vertex
    fx = zeros(size(m{q}.points,1),ns);
    tx = zeros(size(m{q}.points,1),1);
    for j = 1:ncp
       d=sqrt(sum((m{q}.points - ones(size(m{q}.points,1),1) * xcp(j,:)).^2,2)) / rcp(j);
       phi=exp(-d .^ 2);
       fx = fx + phi * fcp(j,:);
       tx = tx + phi * tcp(j);
    end

    % Add some random noise for good measure
    fx = fx + randn(size(fx));

    m1{q} = vtk_add_point_data(m{q}, 'Y', fx);
    m1{q} = vtk_add_point_data(m1{q}, 'truth', tx);
    
    % Now for the cell data
    if q == 1
        ic = cell2mat(m{q}.cells.polygons);
        xc = (m{q}.points(ic(1,:),:) + m{q}.points(ic(2,:),:) + m{q}.points(ic(3,:),:))/3;
    else
        ic = cell2mat(m{1}.cells.cells);
        xc = (m{q}.points(ic(1,:),:) + m{q}.points(ic(2,:),:) + ...
            m{q}.points(ic(3,:),:)+ m{q}.points(ic(4,:),:))/4;
    end
    
    fc = zeros(size(ic,2), ns);
    tc = zeros(size(ic,2),1);
    for j = 1:ncp
       d=sqrt(sum((xc - ones(size(xc,1),1) * xcp(j,:)).^2,2)) / rcp(j);
       phi=exp(-d .^ 2);
       fc = fc + phi * fcp(j,:);
       tc = tc + phi * tcp(j);
    end

    % Add some random noise for good measure
    fc = fc + randn(size(fc));

    m1{q} = vtk_add_cell_data(m1{q}, 'Y', fc);
    m1{q} = vtk_add_cell_data(m1{q}, 'truth', tc);
    
end

%% Save the mesh with fx in it
vtk_polydata_write('test_meshglm_surface.vtk', m1{1});
vtk_mesh_write('test_meshglm_volume.vtk', m1{2});
save('test_meshglm_surface_design.txt', 'X', '-ascii');
save('test_meshglm_surface_contrast.txt', 'con', '-ascii');
