function [xsamp ysamp] = sampcurv(x, y, N, delta)
% SAMPCURV Sample contour with curvature-dependent density
%   [xsam, ysam] = sampcurv(x, y, N, delta) samples the contour (x,y)
%   with N samples. Delta is a parameter that determines how much the 
%   curvature influences sampling density. Delta > 0. Small values mean
%   more influence. Sampling rate is proportional to 1/(curv + delta). 
%   Recommended values are 0.02 - 0.1 or so, but depends on application

t=linspace(0,1,length(x));


dx=gradient(x,t); 
dy=gradient(y,t); 
ddx=gradient(dx,t); 
ddy=gradient(dy,t);
ds=sqrt(dx.*dx+dy.*dy);
curv = (-ddx .* dy + ddy .* dx) ./ (ds.^3);

pp = interp1(t, ds .* (abs(curv) + delta) ,'cubic','pp');

[tode,vode]=ode45(@(t,u)(ppval(pp,t)), t, 0, odeset('RelTol',1e-6));
L = vode(end);

u = linspace(0,L,length(t));
[tode2,vode2] = ode45(@(t,y)(1/ppval(pp,y)), u, 0, odeset('RelTol',1e-6));


tsamp=interp1(t, vode2, linspace(0,1,N+1));
tsamp=tsamp(1:end-1);
xsamp=interp1(t,x,tsamp);
ysamp=interp1(t,y,tsamp);

