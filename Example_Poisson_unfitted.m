clearvars;
close all;
%
% Load the point sets.
%
data = load('data/nodes_unfitted_case.mat');

X.nds = data.X; % X point set. Points enforced in a box larger than the domain Omega.

Y.nds = data.Y; % Y point set. Points that conform to the domain Omega.
Y.idx_dirichlet = data.idx_bc_Y; % Boundary condition indices.
Y.idx_in = data.idx_in_Y; % Interior nodes indices.

M = size(Y.nds,1);

%
% Parameters.
%
bf = Basis.basisF('PHS', '2d'); % Choose basis functions.
p = 3; % PHS power (r^p).
polydeg = 3; % Augmented polynomial degree.
n = 2*nchoosek(polydeg+2,2); % Stencil size.

%
% Remove X points which are too far away from the boundary of Omega.
%
idx_tmp = knnsearch(X.nds,Y.nds,'k',ceil(0.5*n));

X.nds = X.nds(unique(knnsearch(X.nds,Y.nds,'k',ceil(0.5*n))), :);
N = size(X.nds,1);

%
% Compute evaluation/differentiation matrices.
%
[E, Dx, Dy, Dxx, Dyy, Dxy, stencils] = Matrices.generate_2d(X.nds,Y.nds,p,n,polydeg,bf,[]);

%
% Solve a Poisson equation with Dirichlet+Neumann BC data.
%

% Decide for an exact solution.
u_exact = @(x,y) sin(2*pi*x.*y);

% The corresponding right-hand-sides.
f2 = @(x,y) x.^2.*pi.^2.*sin(x.*y.*pi.*2.0).*-4.0-y.^2.*pi.^2.*sin(x.*y.*pi.*2.0).*4.0; % Interior RHS.
f0 = @(x,y) u_exact(x,y); % Dirichlet RHS.

% Assemble the PDE operator.
D = spalloc(M,N,n);
D(Y.idx_in,:) = Dxx(Y.idx_in,:) + Dyy(Y.idx_in,:);
D(Y.idx_dirichlet,:) = E(Y.idx_dirichlet,:);

% Assemble the RHS.
f = zeros(M,1);
f(Y.idx_in) = f2(Y.nds(Y.idx_in,1), Y.nds(Y.idx_in,2)); % Interior.
f(Y.idx_dirichlet) = f0(Y.nds(Y.idx_dirichlet,1), Y.nds(Y.idx_dirichlet,2)); % Dirichlet.

% Scale the PDE operator and the RHS.
[~,dist] = knnsearch(X.nds,X.nds,'k',2); % An approximate distance between interpolation (X) points.
h = mean(dist(:,2));

M0 = length(Y.idx_dirichlet);
M2 = length(Y.idx_in);

D(Y.idx_in,:) = 1/sqrt(M2) * D(Y.idx_in,:);
D(Y.idx_dirichlet,:) = 1/h * 1/sqrt(M0) * D(Y.idx_dirichlet,:);

f(Y.idx_in) = 1/sqrt(M2) * f(Y.idx_in);
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