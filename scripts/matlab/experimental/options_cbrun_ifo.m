%% Set MVU options

% Regularization parameter for MVU (most crucial for good output)
opt.mvu.eta = 1e-02;

% Number of laplacian eigenvectors for MVU (don't mess with it)
opt.mvu.leigsdim = 10;

% Max MVU iterations (don't mess with it)
opt.mvu.maxiter = 2000;

%% Set polynomial fitting options

% Order of polynomial fitting (3 to 5 is reasonable)
opt.fit.order = 5;

%% Set options used to compute 2D contour 

% Size of the 2D raster for embedding in pixels (200 is reasonable)
opt.image.size = 200;

% Size of the closing operator (used to convert scattered points to region)
opt.image.r_closing = 8;

% Size of the dilation operator applied after closing (to avoid sharp
% corners)
opt.image.r_post_dilation = 6;

% Number of iterations of laplacian smoothing applied to the contour
opt.image.smooth_iter = 600;

% Number of boundary samples
opt.contour.n_samples = 90;

% Curvature sensitivity during sampling. Smaller values make sampling
% more proportional to the curvature
opt.contour.curv_sensitivity = 1.2;

%% Set options for the biharmonic cm-rep type

% Whether an inner ring of vertices is created (for biharmonic only)
opt.biharm.inner = 1;

% Initial radius along the cm-rep boundary
opt.biharm.radius = 1;

% Initial value of rho on cm-rep interior
opt.biharm.rho = -0.01;


%% Set options for initial triangulation of contour (unimportant ?)

% Quality of the triangles desired (default 32)
opt.triangle.quality = 32;

% Maximum area for the triangle
opt.triangle.max_area = 60;

% Extra options to triangle
opt.triangle.extra_opts = '';


%% Set the options for AFRONT triangulation of 3D surface

% Smaller values result in more triangles. Default = 0.5
opt.afront.rho = 0.3;