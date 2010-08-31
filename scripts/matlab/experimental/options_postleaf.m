%% Set MVU options

% Regularization parameter for MVU (most crucial for good output)
opt.mvu.eta = 1e-03;

% Number of laplacian eigenvectors for MVU (don't mess with it)
opt.mvu.leigsdim = 10;

% Max MVU iterations (don't mess with it)
opt.mvu.maxiter = 4000;

%% Set polynomial fitting options

% Order of polynomial fitting (3 to 5 is reasonable)
opt.fit.order = 9;

%% Set options used to compute 2D contour 

% Size of the 2D raster for embedding in pixels (200 is reasonable)
opt.image.size = 200;

% Size of the closing operator (used to convert scattered points to region)
opt.image.r_closing = 16;

% Size of the dilation operator applied after closing (to avoid sharp
% corners)
opt.image.r_post_dilation = 0;

% Number of iterations of laplacian smoothing applied to the contour
opt.image.smooth_iter = 200;

% Number of boundary samples
opt.contour.n_samples = 40;

% Curvature sensitivity during sampling. Smaller values make sampling
% more proportional to the curvature (0.02 is a good value)
opt.contour.curv_sensitivity = 0.04;

%% Set options for the cm-rep type

% Whether an inner ring of vertices is created (for biharmonic only)
opt.biharm.inner = 1;

% How far the inner vertices are from parent vertices, fraction of edge
% length (default is 2, i.e. half the edge length)
opt.biharm.disp_factor = 2;

% How the angle is divided when inserting new vertices (default is to
% divide the angle into 3 sectors)
opt.biharm.angle_factor = 3;

%% Set options for triangulation of contour

% Quality of the triangles desired (default 32)
opt.triangle.quality = 32;

% Maximum area for the triangle
opt.triangle.max_area = 60;

% Extra options to triangle
opt.triangle.extra_opts = '';