clearvars;
close all;
%
% Load interpolation points.
%
X = load('Example_Poisson_fitted_nds_N=2000.mat');
N = size(X.nds,1);
%
% Load evaluation points.
%
Y = load('Example_Poisson_fitted_nds_N=6000.mat');
M = size(Y.nds,1);

%
% Parameters.
%
bf = Basis.basisF('PHS', '2d'); % Choose basis functions.
p = 3; % PHS power (r^p).
polydeg = 3; % Augmented polynomial degree.
n = 2*nchoosek(polydeg+2,2); % Stencil size.

%
% Map certain Y points to X points.
%
idx = knnsearch(Y.nds, X.nds, 'k', 1);
Y.nds(idx,:) = X.nds;

%
% Compute the nearest neighbors indeces.
%
indeces.idx_X = knnsearch(X.nds,X.nds,'k',n);
indeces.idx_Y_X = knnsearch(X.nds,Y.nds,'k',1);

%
% Compute evaluation/differentiation matrices.
%
[E, Dx, Dy, Dxx, Dyy, Dxy, stencils] = Matrices.generate_2d(X.nds,Y.nds,p,n,polydeg,bf,indeces,[]);

%
% Solve a Poisson equation with Dirichlet+Neumann BC data.
%

% Decide for an exact solution.
u_exact = @(x,y) sin(2*pi*x.*y);

% The corresponding right-hand-sides.
f2 = @(x,y) x.^2.*pi.^2.*sin(x.*y.*pi.*2.0).*-4.0-y.^2.*pi.^2.*sin(x.*y.*pi.*2.0).*4.0; % Interior RHS.
f1 = @(n1,n2,x,y) n2.*x.*pi.*cos(x.*y.*pi.*2.0).*2.0+n1.*y.*pi.*cos(x.*y.*pi.*2.0).*2.0; % Neumann RHS.
f0 = @(x,y) u_exact(x,y); % Dirichlet RHS.

% Assemble the PDE operator.
D = spalloc(M,N,n);
D(Y.idx_in,:) = Dxx(Y.idx_in,:) + Dyy(Y.idx_in,:);
D(Y.idx_neumann,:) = Y.nrmls(:,1).*Dx(Y.idx_neumann,:) + Y.nrmls(:,2).*Dy(Y.idx_neumann,:);
D(Y.idx_dirichlet,:) = E(Y.idx_dirichlet,:);

% Assemble the RHS.
f = zeros(M,1);
f(Y.idx_in) = f2(Y.nds(Y.idx_in,1), Y.nds(Y.idx_in,2)); % Interior.
f(Y.idx_neumann) = f1(Y.nrmls(:,1), Y.nrmls(:,2), Y.nds(Y.idx_neumann,1), Y.nds(Y.idx_neumann,2)); % Neumann.
f(Y.idx_dirichlet) = f0(Y.nds(Y.idx_dirichlet,1), Y.nds(Y.idx_dirichlet,2)); % Dirichlet.

% Scale the PDE operator and the RHS.
[~,dist] = knnsearch(X.nds,X.nds,'k',2); % An approximate distance between interpolation (X) points.
h = mean(dist(:,2));

M0 = length(Y.idx_dirichlet);
M1 = length(Y.idx_neumann);
M2 = length(Y.idx_in);

D(Y.idx_in,:) = 1/sqrt(M2) * D(Y.idx_in,:);
D(Y.idx_neumann,:) = 1/sqrt(M1) * D(Y.idx_neumann,:);
D(Y.idx_dirichlet,:) = 1/h * 1/sqrt(M0) * D(Y.idx_dirichlet,:);

f(Y.idx_in) = 1/sqrt(M2) * f(Y.idx_in);
f(Y.idx_neumann) = 1/sqrt(M1) * f(Y.idx_neumann);
f(Y.idx_dirichlet) = 1/h * 1/sqrt(M0) * f(Y.idx_dirichlet);

% Solve.
u_Y = E*(D\f);

% Compute error.
err = norm(u_Y - u_exact(Y.nds(:,1),Y.nds(:,2)),2)/norm(u_exact(Y.nds(:,1), Y.nds(:,2)),2)

% Visualize.
figure;
scatter(Y.nds(:,1), Y.nds(:,2), [], u_Y, 'filled'); 
axis equal;
colorbar;