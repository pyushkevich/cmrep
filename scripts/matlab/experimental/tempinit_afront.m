%% READ OPTIONS FIRST
% open file options_XXXX
skel_file = '/Users/pauly/tempinit/posterior_skel_qc_clean.vtk';

%% Load the skeleton from VTK mesh
m=vtk_polydata_read(skel_file);
X=m.points;
R=vtk_get_point_data(m,'Radius');

% Compute the adjacency matrix
adj = vtk_adjacency_matrix(m);

% Get all directed edges
[ei ej]=ind2sub(size(adj), find(adj));

% Plot the points in 3D
clf; scatter3(X(:,1), X(:,2), X(:,3), [], R, '.'); 
axis equal; axis vis3d; colormap jet;

%%  Generate a flat representation of the data

% Compute the squared distance matrix
DD=sparse(ei, ej, sum((m.points(ei,:) - m.points(ej,:)).^2,2));

% Call FAST-MVU on the scattered points
[ymvu detail]=fastmvu(DD, 2, ...
    'leigsdim', opt.mvu.leigsdim,...
    'eta', opt.mvu.eta, 'maxiter', opt.mvu.maxiter);

% Normalize u,v to range 0,1
uv=ymvu';
uv=(uv-ones(length(uv),1)*min(uv))./(ones(length(uv),1)*(max(uv)-min(uv)));
u=uv(:,1); v=uv(:,2);

% Plot the relative distortion in edge length
DDuv=sparse(ei, ej, sum((ymvu(:,ei) - ymvu(:,ej)).^2,1));
scatter(0.5 * (u(ei) + u(ej)), 0.5 * (v(ei) + v(ej)),[],...
    (nonzeros(DDuv) - nonzeros(DD)) ./ nonzeros(DD),'.');
axis equal; axis vis3d; colorbar;

%clf; scatter(u, v);

%% Fit polynomial model to the points
order = opt.fit.order;

% Get polynomial coefficients
beta = tempinit_polyfit(order,[u v], X);
beta_R = tempinit_polyfit(order,[u v], R);

% Visualize the best fit
Xfit = tempinit_polyinterp(order,  beta, [u v]);

% Plot the original data vs. fitted data
scatter3(X(:,1),X(:,2),X(:,3),'bo'); hold on;
scatter3(Xfit(:,1),Xfit(:,2),Xfit(:,3),'r.'); hold off; axis image;

%% Convert representation to an image
np = opt.image.size;
iuv = 1+floor(uv * (np-1));

%I=zeros(np,np);
%I(sub2ind(size(I), iuv(:,1), iuv(:,2)))=1;

I = zeros(np,np);
for i = 1:length(m.cells.polygons)
    I = max(I,poly2mask(...
        iuv(m.cells.polygons{i},1), ...
        iuv(m.cells.polygons{i},2), np, np)');
end

%tri = cell2mat(m.cells.polygons)';
%I=poly2mask(iuv(tri,1),iuv(tri,2),np,np)';

% Pad the image
npad = max(opt.image.r_closing,opt.image.r_post_dilation) + 4;
Ipad = zeros(size(I)+2*npad);
Ipad(npad+1:npad+np, npad+1:npad+np) = I;
imagesc(Ipad'); axis image; colormap gray;

%% Apply morphological closing and smoothing
Iclose = imclose(Ipad, strel('disk',opt.image.r_closing));
Iclose = imdilate(Iclose, strel('disk',opt.image.r_post_dilation));

imagesc(Ipad'); axis image; colormap gray; hold on;
contour(Iclose', [0.5 0.5], 'r'); hold off;

%% Extract and smooth contour
C = contourc(Iclose', [0.5, 0.5]);
ncp = C(2,1);
cuv = C(:,2:ncp+1);

cuv = cuv(:,end:-1:1);

% Smooth the contour using laplace smoothing
for i = 1:opt.image.smooth_iter
    cuv=0.5 * (circshift(cuv,[0 1]) + circshift(cuv,[0 -1]));
end

imagesc(Ipad'); axis image; colormap gray; hold on;
line(cuv(1,:), cuv(2,:));

%% Sample the contour uniformly (in 3 space)
nsam = opt.contour.n_samples;
Xuv=tempinit_polyinterp(order, beta, (cuv' - npad) / np);
[uv_sam X_sam] = tempinit_sampcurv(cuv', Xuv, nsam, ...
    opt.contour.curv_sensitivity, 0);

clf; 
scatter3(X_sam(:,1), X_sam(:,2), X_sam(:,3),'.');
axis equal; axis vis3d;

%% Perform initial 2D triangulation

% Triangulate the sample points
[T N] = tempinit_triangle(uv_sam, [], opt.triangle.quality, opt.triangle.max_area);

% Map them to 3-space
Xn = tempinit_polyinterp(order, beta, (N - npad) / np);
Rn = tempinit_polyinterp(order, beta_R, (N - npad) / np);

clf;
imagesc(Ipad'); axis image; colormap jet; hold on;
trisurf(T,N(:,1),N(:,2),Rn,'EdgeColor','white','FaceColor','interp');

%% Run AFRONT

% Export to an .off format mesh (what seems to work)
f = fopen('initmesh.off', 'w+');
fprintf(f,'OFF\n%i %i 0\n',size(Xn,1), size(T,1));
fprintf(f,'%f %f %f\n', Xn');
fprintf(f,'3 %i %i %i\n', T'-1);
fclose(f);

% Run AFRONT
rho = opt.afront.rho;
cmd = sprintf('afront -nogui initmesh.off -rho %i -tri_mesh 1', rho);
system(cmd);

% Read output file
Taf = []; Xaf = [];
f = fopen('outmesh.m','rt');
while ~ feof(f)
    ln = fgets(f);
    if ln(1) == 'V'
        data = sscanf(ln,'Vertex %i %f %f %f\n');
        Xaf(data(1),:) = data(2:4);
    else
        data = sscanf(ln,'Face %i %i %i %i\n');
        Taf(data(1),:) = data(2:4);
    end
end
fclose(f);

clf;
trisurf(Taf, Xaf(:,1), Xaf(:,2), Xaf(:,3), 'EdgeColor','black','FaceColor','interp');
axis equal; axis vis3d;


%% Split triangles at the edge for biharmonic models

tr = TriRep(Taf, Xaf);

% Sparse matrix representing what edges have been split
eall = edges(tr);
esplit = sparse(eall(:,1), eall(:,2), -1);

% Get the list of all triangles on the boundary
ebnd = freeBoundary(tr);
fbnd = edgeAttachments(tr,ebnd);

% Counter for inserted vertices
vlast = size(Xaf,1);

% Storage for new vertices
Xsp = nan(vlast+size(ebnd,1)*3,3);
Xsp(1:vlast,:) = Xaf;

% Storage for new triangles
Tsp = [];
tlast = 0;

% Current triangles marked for deletion
tdel = zeros(size(Taf,1),1);

% Loop over boundary edges
for i = 1:size(ebnd,1)
   
    eisrt = ebnd(i,:);
    ti = cell2mat(edgeAttachments(tr,eisrt));
    
    vopp = setdiff(Taf(ti,:), eisrt);
    iopp = find(Taf(ti,:)==vopp);
    ei = [Taf(ti,mod(iopp,3)+1) Taf(ti,mod(iopp+1,3)+1)]; 
    
    % Position of edge, opposite vertex in triangle
    idx=[find(Taf(ti,:)==ei(1)), find(Taf(ti,:)==ei(2)), find(Taf(ti,:)==vopp)]; 
    
    % Lengths of opposing edges in triangle
    l = sqrt(sum((Xaf(circshift(Taf(ti,:),[0 1]),:) - ...
        Xaf(circshift(Taf(ti,:),[0 -1]),:)).^2,2));
    
    % Split existing edges 
    vcut = [0 0];
    for j = 1:2
        vcut(j) = esplit(min(vopp,ei(j)), max(vopp,ei(j)));
        if vcut(j) < 0
            vcut(j) = vlast + 1;
            vlast = vcut(j);
            esplit(min(vopp,ei(j)), max(vopp,ei(j))) = vcut(j);
            
            xbary=[0,0,0];
            xbary(idx(3))=l(idx(3));
            xbary(idx(j))=l(idx(j));
            Xsp(vcut(j),:) = baryToCart(tr,ti,xbary./sum(xbary));
        end
    end
    
    % Add center vertex
    vlast = vlast + 1;
    Xsp(vlast,:) = baryToCart(tr,ti,l'./sum(l));
    
    % Add new triangles
    Tsp(tlast+1,:) = [vlast vcut(2) vcut(1)];
    Tsp(tlast+2,:) = [vlast ei(2) vcut(2)];
    Tsp(tlast+3,:) = [vlast ei(1) ei(2)];
    Tsp(tlast+4,:) = [vlast vcut(1) ei(1)];
    Tsp(tlast+5,:) = [vopp vcut(1) vcut(2)];
    tlast = tlast + 5;
    
    % Mark old triangle for deletion
    tdel(ti) = 1;
end 

% Loop over all other triangles
for i = 1:size(Taf,1)
    % Ignore processed triangles
    if tdel(i), continue; end
    
    % Check if any of the edges were cut
    tri = Taf(i,:);
    vcut = [...
        esplit(min(tri(2),tri(3)), max(tri(2),tri(3))),...
        esplit(min(tri(1),tri(3)), max(tri(1),tri(3))),...
        esplit(min(tri(2),tri(1)), max(tri(2),tri(1)))];
    
    % Compute angles
    dx1 = Xaf(circshift(tri,[0 -1]),:) - Xaf(tri,:);
    dx2 = Xaf(circshift(tri,[0  1]),:) - Xaf(tri,:);
    cosj = dot(dx1',dx2') ./ sqrt(dot(dx1',dx1') .* dot(dx2',dx2'));
    
    % Add new triangles
    ncut = sum(vcut > 0);
    if ncut == 0
        continue;
    elseif ncut == 1
        q = find(vcut > 0);
        q1 = mod(q,3)+1; q2 = mod(q+1,3)+1; 
        Tsp(tlast + 1,:) = [tri(q), tri(q1), vcut(q)];
        Tsp(tlast + 2,:) = [vcut(q), tri(q2), tri(q)];
        tlast = tlast + 2;
    elseif ncut == 2
        q = find(vcut < 0);
        q1 = mod(q,3)+1; q2 = mod(q+1,3)+1; 
        Tsp(tlast + 1,:) = [tri(q), vcut(q2), vcut(q1)];
        if cosj(q1) < cosj(q2)
            Tsp(tlast + 2,:) = [tri(q1) vcut(q1) vcut(q2)];
            Tsp(tlast + 3,:) = [tri(q1) tri(q2) vcut(q1)];
        else
            Tsp(tlast + 2,:) = [tri(q2) vcut(q1) vcut(q2)];
            Tsp(tlast + 3,:) = [tri(q2) vcut(q2) tri(q1)];
        end
        tlast = tlast + 3;
    else
        Tsp(tlast + 1,:) = [tri(q) vcut(q2) vcut(q1)];
        Tsp(tlast + 2,:) = [tri(q1) vcut(q) vcut(q2)];
        Tsp(tlast + 3,:) = [tri(q2) vcut(q1) vcut(q)];
        Tsp(tlast + 4,:) = [vcut(q) vcut(q1) vcut(q2)];
        tlast = tlast + 1;
    end
    tdel(i) = 1;
end

Xsp=Xsp(all(isfinite(Xsp),2),:);
Tsp=[Taf(tdel==0,:); Tsp];

clf;
trisurf(Tsp, Xsp(:,1), Xsp(:,2), Xsp(:,3), 'EdgeColor','black','FaceColor','interp');
axis equal; axis vis3d;

%% Export the model using VTK
mnew.hdr=m.hdr;
mnew.points=Xsp;
mnew.cells.polygons=num2cell(Tsp',1);
if(opt.biharm.inner ~= 1)
    mnew=vtk_add_point_data(mnew,'Radius',R2 / 2,1);
end
vtk_polydata_write('medialtemplate.vtk', mnew);

% Write the cmrep file - you may want to edit this for different models
f = fopen('medialtemplate.cmrep','wt');
fprintf(f,'Grid.Type = LoopSubdivision\n');
fprintf(f,'Grid.Model.SolverType = PDE\n');
fprintf(f,'Grid.Model.Atom.SubdivisionLevel = 0\n');
fprintf(f,'Grid.Model.Coefficient.FileName = medialtemplate.vtk\n');
fprintf(f,'Grid.Model.Coefficient.FileType = VTK\n');

if(opt.biharm.inner == 1)
    fprintf(f,'Grid.Model.Coefficient.ConstantRho = %f\n',...
        opt.biharm.rho);
    fprintf(f,'Grid.Model.Coefficient.ConstantRadius.Boundary = %f\n',...
        opt.biharm.radius);
    fprintf(f,'Grid.Model.SolverType = PDE\n');
else
    fprintf(f,'Grid.Model.SolverType = BruteForce\n');
end
    
fclose(f);

