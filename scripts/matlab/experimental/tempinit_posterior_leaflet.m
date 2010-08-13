%% Load the skeleton from VTK mesh
m=vtk_polydata_read('posterior_skel.vtk');
rad=vtk_get_point_data(m,'Radius');

% Randomly sample ~1000 points
qperm=randperm(size(m.points,1));
qsamp=qperm(1:400);

% Extract the X and R of these points
X=m.points(qsamp,:);
R=rad(qsamp);

scatter3(X(:,1), X(:,2), X(:,3), [], R); axis vis3d;

%%  Generate a flat representation of the data
% X=load('-ascii', 'skel_xyz.mat'); 

% Call MVU on the scattered points
uv = mvuIncXL(X', 6, 400, 2, struct('factor',0.99999))';

% Normalize u,v to range 0,1
uv=(uv-ones(length(uv),1)*min(uv))./(ones(length(uv),1)*(max(uv)-min(uv)));
u=uv(:,1); v=uv(:,2);
clf; scatter(u, v);

%% Fit polynomial model to the points

% Fit fourth-order polynomial model X = F * beta + eps
order=3;
F=zeros(length(u), order*(order+1)/2);
for j=0:order
    for k=0:j
        F(:,1+j*(j+1)/2+k)=(u.^k) .* (v.^(j-k));
    end
end

beta=F\X;
Xfit = F*beta;

beta_R=F\R;

% Plot the original data vs. fitted data
scatter3(X(:,1),X(:,2),X(:,3),'bo'); hold on;
scatter3(Xfit(:,1),Xfit(:,2),Xfit(:,3),'r.'); hold off; axis vis3d;

%% Convert representation to an image
np = 200;
iuv = 1+floor(uv * (np-1));
I=zeros(np,np);
I(sub2ind(size(I), iuv(:,1), iuv(:,2)))=1;
imagesc(I);

% Pad the image
npad=20;
Ipad = zeros(size(I)+2*npad);
Ipad(npad+1:npad+np, npad+1:npad+np) = I;
imagesc(Ipad'); axis image; colormap gray;

%% Apply morphological closing and smoothing
Iclose = imclose(Ipad, strel('disk',16));
Iclose = imdilate(Iclose, strel('disk',12));

Ismooth = imfilter(Iclose, fspecial('gaussian', 96, 12.0));
imagesc(Ipad'); axis image; colormap gray; hold on;
contour(Iclose', [0.5 0.5], 'r'); hold off;
%contour(Ismooth', [0.5 0.5], 'r'); hold off;

%% Extract and smooth contour
C = contourc(Iclose', [0.5, 0.5]);
ncp = C(2,1);
cuv = C(:,2:ncp+1);

% Smooth the contour using laplace smoothing
for i = 1:800
    cuv=0.5 * (circshift(cuv,[0 1]) + circshift(cuv,[0 -1]));
end

imagesc(Ipad'); axis image; colormap gray; hold on;
line(cuv(1,:), cuv(2,:));

%% Sample the contour

nsam = 91;
[xsam ysam] = sampcurv(cuv(1,:), cuv(2,:), nsam, 0.04);
csam=[xsam; ysam];

imagesc(Ipad'); axis image; colormap gray; hold on;
line(csam(1,:), csam(2,:));
scatter(csam(1,:), csam(2,:), 'b.');





%% For biharmonic models, construct an inner contour
biharm=1;
if(biharm)
    % We want to generate two new samples for every boundary vertex. There
    % must be many ways of doing this
    A = csam; Al = circshift(A, [0 1]); Ar = circshift(A, [0 -1]);
    vl = (Al - A) ./ ([1; 1] * sqrt(diag((Al-A)' * (Al-A)))');
    vr = (Ar - A) ./ ([1; 1] * sqrt(diag((Ar-A)' * (Ar-A)))');
    theta = acos(diag(vl' * vr)) / 3;
    rl = diag((Al - A)' * vl) / 2; 
    rr = diag((Ar - A)' * vr) / 2;

    ivl = A + ([1;1] * rl') .* ...
       [cos(theta)' .* vl(1,:) - sin(theta)' .* vl(2,:);
       sin(theta)' .* vl(1,:) + cos(theta)' .* vl(2,:)];

    ivr = A + ([1;1] * rr') .* ...
       [cos(theta)' .* vr(1,:) + sin(theta)' .* vr(2,:);
       -sin(theta)' .* vr(1,:) + cos(theta)' .* vr(2,:)];

    imagesc(Ipad'); axis image; colormap gray;
    line(csam(1,:), csam(2,:))
    scatter(csam(1,:), csam(2,:), 'b.');
    scatter(ivl(1,:), ivl(2,:), 'r.');
    scatter(ivr(1,:), ivr(2,:), 'g.');

    nextra = nsam * 2;
else
    nextra = 0;
end

%% Run triangle

% Generate poly file
f = fopen('contour.poly', 'w+');
fprintf(f,'%i 2 0 0\n', nsam + nextra);
fprintf(f,'%i %f %f\n', [1:nsam; csam]);
fprintf(f,'%i %f %f\n', [nsam+1:nsam+nextra; horzcat(ivl, ivr)]);
fprintf(f,'%i 0\n', nsam);
fprintf(f,'%i %i %i\n', [1:nsam; 1:nsam; 1+mod([1:nsam],nsam)]);
fprintf(f,'0\n');
fclose(f);

% Run triangle
system('triangle -Y -pq32 -a120 contour.poly');

%% Read its output

% Read element file
f = fopen('contour.1.ele','rt');
head = fscanf(f, '%i',3);
body = fscanf(f, '%i',[4 head(1)])';
T=body(:,2:4);
fclose(f);

% Read node file
f = fopen('contour.1.node','rt');
head = fscanf(f, '%i',4);
body = fscanf(f, '%f',[4 head(1)])';
N=body(:,2:3);
fclose(f);

imagesc(Ipad'); axis image; colormap gray; hold on;
trimesh(T,N(:,1),N(:,2),'Color','white');
if(biharm)
    scatter(csam(1,:), csam(2,:), 'b.');
    scatter(ivl(1,:), ivl(2,:), 'r.');
    scatter(ivr(1,:), ivr(2,:), 'g.')
end
    


%% Map the vertices to the 3-space

% Map uv back to 0-1 range
Nnorm=(N-npad)/np;
un=Nnorm(:,1); vn=Nnorm(:,2);

% Interpolate the x,y,z of the nodes
Fn=zeros(length(un), order*(order+1)/2);
for j=0:order
    for k=0:j
        Fn(:,1+j*(j+1)/2+k)=(un.^k) .* (vn.^(j-k));
    end
end

Xn = Fn * beta;
Rn = Fn * beta_R;

% Plot the new nodes
clf;
scatter3(X(:,1),X(:,2),X(:,3),'bo'); hold on;
trisurf(T, Xn(:,1), Xn(:,2), Xn(:,3),Rn); hold off; colormap jet;
axis image;

%% Export the model using VTK
mnew.hdr=m.hdr;
mnew.points=Xn;
mnew.cells.polygons=num2cell(T',1);
%mnew=vtk_add_point_data(mnew,'Radius',Rn / 2,1);
vtk_polydata_write('medialtemplate.vtk', mnew);

%% Write the cmrep file - you may want to edit this for different models
f = fopen('medialtemplate.cmrep','wt');
fprintf(f,'Grid.Type = LoopSubdivision\n');
fprintf(f,'Grid.Model.SolverType = PDE\n');
fprintf(f,'Grid.Model.Atom.SubdivisionLevel = 0\n');
fprintf(f,'Grid.Model.Coefficient.FileName = medialtemplate.vtk\n');
fprintf(f,'Grid.Model.Coefficient.FileType = VTK\n');  
fprintf(f,'Grid.Model.Coefficient.ConstantRho = -0.01\n');
fprintf(f,'Grid.Model.Coefficient.ConstantRadius = 0.4\n');
fclose(f);


