% (c) Igor Tominec, Uppsala University
% This function generates the RBF-FD matrices in two dimensions. RBFs are augmented with the polynomial basis.
%
% Input:
% X ... Nx1 list of interpolation points.
% Y ... Mx1 list of evaluation points. (If Y=X, then the matrices will be square)
% epOrP ... exponent on PHS or shape parameter if using GA or MQ basis functions.
% n ... stencil size.
% polydeg ... polynomial degree.
% basis ... the basis functions.
% indeces ... a struct with indeces of stencil neighbors (idx_X) and indeces of closest stencils to every evaluation points (idx_Y_X)
% stencils ... a struct of already computed stencils (can be used for interpolation purposes). Input [] when no precomputed stencils are available.
%
% Output:
% E ... MxN evaluation matrix.
% Dx ... MxN differentiation matrix d/dx
% Dy ... MxN differentiation matrix d/dy
% Dxx ... MxN differentiation matrix d^2/dx^2
% Dxx ... MxN differentiation matrix d^2/dy^2
% Dxy ... MxN differentiation matrix d^2/(dx*dy)

function [Interp, Dx, Dxx, stencils] = generate_1d(X, Y, epOrP, n, polydeg, basis, stencils, varargin)

import Basis.*
import Basis.Matrices.Tools.*


%
% Preallocate the dense matrices where every stencil is put in.
%

dim = size(X,2);


%
% Find the neighbours
%
idx_X = knnsearch(X,X,'k',n);
idx_Y_X = knnsearch(X,Y,'k',1);

idx_centers = unique(idx_Y_X);
idx_centers_local = 1:length(idx_centers);
idx_centers_global = sparse(size(X,1),1);
idx_centers_global(idx_centers) = idx_centers_local;


N_dom = size(X,1);
N_centers = size(idx_centers,1); % X(idx_domain) could be interior, or some of the boundaries.

%
% Precompute the monomial matrix.
%
monMatrix = polyMat_monmatrix(polydeg, dim);

%
% Precompute the differentiated monomial terms. They are equal for all
% stencils.
%
M_inv = cell(N_centers,1);

Dbf = basis.Dbf;
DDbf = basis.DDbf;
bf = basis.bf;
r = basis.r;

% What is the scaling for the given stencil?
X_scaling = zeros(N_dom,dim);
X_scaling_local = zeros(N_dom,dim);


if (isempty(stencils))
    for i = 1:N_centers
        idx_c = idx_centers(i,1);
        idx = idx_X(idx_c,:);
        
        X_local = X(idx,:);
        
        [X_scaled, X_scaling] = Matrices.Tools.scaleStencil(X_local, X(idx_c,:));
        
        [~, Minv_local] = Matrices.Tools.generateInterpMatrix(X_scaled, basis, monMatrix, epOrP, dim);
        %condM(i) = cond(Minv_local);
        M_inv{i} = Minv_local;
        X_scaling_local(i,:) = X_scaling(1,:);
    end
    X_scaling(idx_centers,:) = X_scaling_local(1:N_centers,:);
    
    
    % Save the stencils
    stencils.M_inv = M_inv;
    stencils.X_scaling = X_scaling;
else
    M_inv = stencils.M_inv;
    X_scaling = stencils.X_scaling;
end
%max(condM)

stencils.idx_centers = idx_centers;

X_scaling(idx_centers,:) = X_scaling_local(1:N_centers,:);

% Precomputation: Evaluate the polynomials and its derivatives.

Y_shift_scale = (Y - X(idx_Y_X,:)) .* X_scaling(idx_Y_X,:);
% Y_shift_scale = (Y) .* X_scaling(idx_Y_X,:);
c = polyMat_eval(Y_shift_scale, monMatrix, dim);
cx = polyMat_diff([1], Y_shift_scale, monMatrix, dim);
cxx = polyMat_diff([2], Y_shift_scale, monMatrix, dim);
n_mon = length(cx(1,:));
idx_st = cumsum(repmat(n+n_mon,3,1)); 

% The matrices containing the weights. Not global yet.
N_dom_Y = size(Y,1);
Interp_loc = zeros(N_dom_Y,n);
Dx_loc = zeros(N_dom_Y,n);
Dxx_loc = zeros(N_dom_Y,n);

% If you use parfor here, it will make things slower! (broadcast stuff...)
for k=1:N_dom_Y

    % Take the interpolation matrix from x in X which is closest to the current y in Y.
    idx_inv = idx_Y_X(k);
    idx_inv_local = idx_centers_global(idx_inv);
    
    % The neighbours around the X-center.
    idx = idx_X(idx_inv,:);

    % The Y-center.
    idx_c = k;
    
    %
    % Shift and scale the stencil nodes to a circle-unit.
    %
    Xscaled = zeros(length(idx),dim);
    Yscaled = zeros(1,dim);
    scale = zeros(dim,1);    
    
    % Shift the stencil nodes so that the center node is the origin.
    for j=1:dim
        Xscaled(:,j) = X(idx,j) - X(idx_inv,j);
        Yscaled(j) = Y(idx_c,j) - X(idx_inv,j);
        
        scale(j) = X_scaling(idx_inv,j);
        
        Xscaled(:,j) = Xscaled(:,j)*scale(j);
        Yscaled(j) = Yscaled(j)*scale(j);
        
        
    end
    
    %
    % Stencil distance matrices. Center point is in the origin (0).
    %
    dx_stencil = Yscaled(1) - Xscaled(:,1);

    % Put the differences to eps to avoid divisions by 0.
    if (dx_stencil(1) == 0)
        dx_stencil(1) = eps;
    end

    r_stencil = r(dx_stencil); % here the distance matrix is actually not scaled.
    
    %
    % 0th derivatives
    %
    b = bf(r_stencil,epOrP); % Here the distance is actually not scaled.        

    stencils.b{k} = [b; c(k,:)'];
    %
    % 1st derivatives
    %    
    bx = Dbf(r_stencil, epOrP, dx_stencil);

    %
    % Second derivatives
    %
    bxx = DDbf(r_stencil, epOrP, dx_stencil);
    %
    % Compute the stencils at once.
    %
    stenc =  M_inv{idx_inv_local} * [bx bxx b; cx(k,:)' cxx(k,:)' c(k,:)'];
    
    Dx_loc(k, :) = scale(1).*stenc(1:idx_st(1)-n_mon);
    Dxx_loc(k, :) = scale(1)^2.*stenc(idx_st(1)+1:idx_st(2)-n_mon);
    Interp_loc(k,:) = stenc(idx_st(2)+1:idx_st(3)-n_mon);
    
end
idx_rows = repmat(1:N_dom_Y, n, 1)';
Dx = sparse(idx_rows, idx_X(idx_Y_X,:), Dx_loc, N_dom_Y, N_dom, N_dom_Y*n);
Dxx = sparse(idx_rows, idx_X(idx_Y_X,:), Dxx_loc, N_dom_Y, N_dom, N_dom_Y*n);
Interp = sparse(idx_rows, idx_X(idx_Y_X,:), Interp_loc, N_dom_Y, N_dom, N_dom_Y*n);


% figure;
% plot(log10(condM), '.');
% title(['log10(cond(M))' ', N=' num2str(size(X,1)) ]);
