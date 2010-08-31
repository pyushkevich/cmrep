function [T N] = tempinit_triangle(uv, uvextra, qual, maxarea, extraopts)
% tempinit_triangle: run TRIANGLE to create a 2D nice mesh
%     [T N] = tempinit_triangle(uv, uvextra, opt) creates a triangulation
%     from contour points uv (n by 2) and optional internal points
%     uvextra. 

nsam = size(uv, 1);
nextra = size(uvextra, 1);

if nargin < 5, extraopts = ''; end;
if nargin < 4, maxarea = 1e100; end;
if nargin < 3, qual = 32; end;

% Generate poly file
f = fopen('contour.poly', 'w+');
fprintf(f,'%i 2 0 0\n', nsam + nextra);
fprintf(f,'%i %f %f\n', [1:nsam; uv']);
if ~isempty(uvextra)
    fprintf(f,'%i %f %f\n', [nsam+1:nsam+nextra; uvextra']);
end
fprintf(f,'%i 0\n', nsam);
fprintf(f,'%i %i %i\n', [1:nsam; 1:nsam; 1+mod([1:nsam],nsam)]);
fprintf(f,'0\n');
fclose(f);

% Run triangle
system(sprintf('triangle -Y -pq%f -a%f %s contour.poly',...
    qual, maxarea, extraopts));

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