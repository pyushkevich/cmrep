function [uv_samp X_samp] = tempinit_sampcurv(uv, X, N, delta, sm_iter)
% SAMPCURV Sample parametrized contour with curvature-dependent density
%   [uv_samp X_samp] = tempinit_sampcurv(uv, X, N, delta, sm_iter) samples 
%   contour X=f(uv) with N samples. 
%   
%   Delta is a parameter that determines how much the 
%   curvature influences sampling density. Delta > 0. Small values mean
%   more influence. Sampling rate is proportional to 1/(curv + delta). 
%   Recommended values are 0.02 - 0.1 or so, but depends on application
%
%   sm_iter is the number of iterations of laplacian smoothing applied to
%   the contour before curvature is computed. For more rough contours,
%   values 100 - 500 are recommended.

n = length(uv);

% Smooth the contour
Y=X;
for i=1:sm_iter, 
    Y=0.5*(circshift(Y,[1,0]) + circshift(Y,[-1,0])); 
end

% Compute the curvature of the smoothed contour
h = 1/n;
dY = (circshift(Y,[-1 0]) - circshift(Y,[1 0])) / (2 * h);
ddY = (circshift(dY,[-1 0]) - circshift(dY,[1 0])) / (2 * h);
ds=sqrt(sum(dY.^2,2));
curv = sqrt(sum(cross(dY, ddY).^2,2)) ./ (ds.^3);

% Compute curvature for every line segment
curvseg = 0.5 * (curv + circshift(curv,-1));

% Compute segment length for every segment
slen = sum((X-circshift(X,[-1 0])).^2,2).^.5;
slenscale = slen .* (curvseg + delta);

% Compute cumulative segment length, and invert
slencum = [0; cumsum(slenscale)];
L=slencum(end);

% Compute inverse via linear interpolation
t = linspace(0,1,n+1);
tsamp = interp1(slencum,t,linspace(0,L,N+1));

% Interpolate uv and X
uv_samp=interp1(t, [uv; uv(1,:)], tsamp);
X_samp=interp1(t,[X; X(1,:)], tsamp);

% Drop the last sample (because it equals the first)
uv_samp = uv_samp(1:end-1,:);
X_samp = X_samp(1:end-1,:);

%pp = interp1(t, ds .* (abs(curv) + delta) ,'cubic','pp');

%[tode,vode]=ode45(@(t,u)(ppval(pp,t)), t, 0, odeset('RelTol',1e-6));
%L = vode(end);

%u = linspace(0,L,length(t));
%[tode2,vode2] = ode45(@(t,y)(1/ppval(pp,y)), u, 0, odeset('RelTol',1e-6));

%tsamp=interp1(t, vode2, linspace(0,1,N+1));
%tsamp=tsamp(1:end-1);


