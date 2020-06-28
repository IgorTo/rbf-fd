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


function [E, Dx, Dy, Dxx, Dyy, Dxy, stencils] = generate_2d(X, Y, epOrP, n, polydeg, basis, stencils, varargin)

import Basis.*
import Basis.Matrices.Tools.*

dim = size(X,2);

%
% Find the neighbours
%
indeces.idx_X = knnsearch(X,X,'k',n);
indeces.idx_Y_X = knnsearch(X,Y,'k',1);

idx_X = indeces.idx_X;
idx_centers = unique(indeces.idx_Y_X);
idx_centers_local = 1:length(idx_centers);
idx_centers_global = sparse(size(X,1),1);
idx_centers_global(idx_centers) = idx_centers_local;

N_dom = size(X,1);
N_centers = length(idx_centers);

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
DDbf_ij = basis.DDbf_ij;
bf = basis.bf;
r = basis.r;

% Allocate storage for stencil scaling.
X_scaling = zeros(N_dom,dim);
X_scaling_local = zeros(N_dom,dim);


if (isempty(stencils))
    for i = 1:N_centers
        idx_c = idx_centers(i);
        idx = idx_X(idx_c,:);
        
        X_local = X(idx,:);
        
        [X_scaled, X_scaling] = Matrices.Tools.scaleStencil(X_local, X(idx_c,:));
        
        [~, Minv_local] = Matrices.Tools.generateInterpMatrix(X_scaled, basis, monMatrix, epOrP, dim);
        
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

%
% Precomputation: Evaluate the polynomials and its derivatives.
%
idx_Y_X = indeces.idx_Y_X;

 % Shift and scale the Y points.
Y_shift_scale = (Y - X(idx_Y_X,:)) .* X_scaling(idx_Y_X,:);

% Compute the polynomial vectors.
c = polyMat_eval(Y_shift_scale, monMatrix, dim);
cx = polyMat_diff([1 0], Y_shift_scale, monMatrix, size(X,2));
cy = polyMat_diff([0 1], Y_shift_scale, monMatrix, size(X,2));
cxx = polyMat_diff([2 0], Y_shift_scale, monMatrix, size(X,2));
cyy = polyMat_diff([0 2], Y_shift_scale, monMatrix, size(X,2));
cxy = polyMat_diff([1 1], Y_shift_scale, monMatrix, size(X,2));
n_mon = length(cx(1,:));
idx_st = cumsum(repmat(n+n_mon,6,1)); 

% The matrices containing the weights. Not global yet.
N_dom_Y = size(Y,1);
E_loc = zeros(N_dom_Y,n);
Dx_loc = zeros(N_dom_Y,n);
Dy_loc = zeros(N_dom_Y,n);
Dxx_loc = zeros(N_dom_Y,n);
Dyy_loc = zeros(N_dom_Y,n);
Dxy_loc = zeros(N_dom_Y,n);

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
    Xscaled = zeros(length(idx),3);
    Yscaled = zeros(1,3);
    
    % Shift the stencil nodes so that the center node is the origin.
    Xscaled(:,1) = X(idx,1) - X(idx_inv,1);
    Xscaled(:,2) = X(idx,2) - X(idx_inv,2);
    
    % Do the same for the Y-center.
    Yscaled(1) = Y(idx_c,1) - X(idx_inv,1);
    Yscaled(2) = Y(idx_c,2) - X(idx_inv,2);

    % Scale the points to [0,1]x[0,1]x...x[0,1]
    scale = zeros(2,1);
    
    scale(1) = X_scaling(idx_inv,1);
    scale(2) = X_scaling(idx_inv,2);

    Xscaled(:,1) = Xscaled(:,1)*scale(1);
    Xscaled(:,2) = Xscaled(:,2)*scale(2);
    
    Yscaled(1) = Yscaled(1)*scale(1);
    Yscaled(2) = Yscaled(2)*scale(2);
    
    %
    % Stencil distance matrices. Center point is in the origin (0).
    %
    dx_stencil = Yscaled(1) - Xscaled(:,1);
    dy_stencil = Yscaled(2) - Xscaled(:,2);

    % Put the differences to eps to avoid divisions by 0.
    if (dx_stencil(1) == 0)
        dx_stencil(1) = eps;
    end
    
    if(dy_stencil(1) == 0)
        dy_stencil(1) = eps;
    end
    
    r_stencil = r(dx_stencil, dy_stencil); % here the distance matrix is actually not scaled.
    
    %
    % 0th derivatives
    %
    b = bf(r_stencil,epOrP); % Here the distance is actually not scaled.        

    %
    % 1st derivatives
    %    
    bx = Dbf(r_stencil, epOrP, dx_stencil);
    by = Dbf(r_stencil, epOrP, dy_stencil);

    %
    % Second derivatives
    %
    bxx = DDbf(r_stencil, epOrP, dx_stencil);
    byy = DDbf(r_stencil, epOrP, dy_stencil);
    bxy = DDbf_ij(r_stencil, epOrP, dx_stencil, dy_stencil);
    
    %
    % Compute the stencils at once.
    %
    stenc =  M_inv{idx_inv_local} * [bx by bxx byy bxy b; cx(k,:)' cy(k,:)' cxx(k,:)' cyy(k,:)' cxy(k,:)' c(k,:)'];
    
    Dx_loc(k, :) = scale(1).*stenc(1:idx_st(1)-n_mon);
    Dy_loc(k, :) = scale(2).*stenc(idx_st(1)+1:idx_st(2)-n_mon);
    Dxx_loc(k, :) = scale(1)^2.*stenc(idx_st(2)+1:idx_st(3)-n_mon);
    Dyy_loc(k, :) = scale(2)^2.*stenc(idx_st(3)+1:idx_st(4)-n_mon);
    Dxy_loc(k, :) = scale(1)*scale(2).*stenc(idx_st(4)+1:idx_st(5)-n_mon);
    E_loc(k,:) = stenc(idx_st(5)+1:idx_st(6)-n_mon);
    
end
idx_rows = repmat(1:N_dom_Y, n, 1)';
Dx = sparse(idx_rows, idx_X(idx_Y_X,:), Dx_loc, N_dom_Y, N_dom, N_dom_Y*n);
Dy = sparse(idx_rows, idx_X(idx_Y_X,:), Dy_loc, N_dom_Y, N_dom, N_dom_Y*n);
Dxx = sparse(idx_rows, idx_X(idx_Y_X,:), Dxx_loc, N_dom_Y, N_dom, N_dom_Y*n);
Dyy = sparse(idx_rows, idx_X(idx_Y_X,:), Dyy_loc, N_dom_Y, N_dom, N_dom_Y*n);
Dxy = sparse(idx_rows, idx_X(idx_Y_X,:), Dxy_loc, N_dom_Y, N_dom, N_dom_Y*n);
E = sparse(idx_rows, idx_X(idx_Y_X,:), E_loc, N_dom_Y, N_dom, N_dom_Y*n);