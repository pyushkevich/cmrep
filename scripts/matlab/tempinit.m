%%  Generate a flat representation of the data
X=load('-ascii', 'skel_xyz.mat');

% Call MVU on the scattered points
uv = mvuIncXL(X', 6, 400, 2, struct('factor',0.99))';

% Normalize u,v to range 0,1
uv=(uv-ones(length(uv),1)*min(uv))./(ones(length(uv),1)*(max(uv)-min(uv)));
u=uv(:,1); v=uv(:,2);
scatter(u, v);

%% Fit polynomial model to the points

% Fit fourth-order polynomial model X = F * beta + eps
order=5;
F=zeros(length(u), order*(order+1)/2);
for j=0:order
    for k=0:j
        F(:,1+j*(j+1)/2+k)=(u.^k) .* (v.^(j-k));
    end
end

beta=F\X;
Xfit = F*beta;

% Plot the original data vs. fitted data
scatter3(X(:,1),X(:,2),X(:,3),'bo'); hold on;
scatter3(Xfit(:,1),Xfit(:,2),Xfit(:,3),'r.'); hold off; axis image;

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
imagesc(Ipad);

%% Apply morphological closing and smoothing
Iclose = imclose(Ipad, strel('disk',16));
Ismooth = imfilter(Iclose, fspecial('gaussian', 24, 4.0));
imagesc(Ipad); axis image; colormap gray; hold on;
contour(Ismooth, [0.5 0.5], 'r'); hold off;

%% Extract and uniformly sample contour
C = contourc(Ismooth', [0.5, 0.5]);
ncp = C(2,1);
cuv = C(:,2:1+ncp);

nsam = 60;
isam = round(linspace(1, ncp, 60));
csam = cuv(:,isam);

scatter(u,v); axis image; colormap gray;
line(csam(1,:), csam(2,:))


%% Run triangle

% Generate poly file
f = fopen('contour.poly', 'wt');
fprintf(f,'%i 2 0 0\n', nsam);
fprintf(f,'%i %f %f\n', [1:nsam; csam]);
fprintf(f,'%i 0\n', nsam);
fprintf(f,'%i %i %i\n', [1:nsam; 1:nsam; 1+mod([1:nsam],nsam)]);
fprintf(f,'0\n');
fclose(f);

% Run triangle
system('triangle -pq32 -a80 contour.poly');

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

imagesc(Ipad); axis image; colormap gray; hold on;
trimesh(T,N(:,1),N(:,2),'Color','red');


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

% Plot the new nodes
clf;
scatter3(X(:,1),X(:,2),X(:,3),'bo'); hold on;
trisurf(T, Xn(:,1), Xn(:,2), Xn(:,3)); hold off;



