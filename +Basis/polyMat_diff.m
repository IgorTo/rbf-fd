function c = polyMat_diff(dX, Xcenter, mon, dim)
% This function computes the vector of the differentiated monomial terms
% evaluated in the center point of the stencil.
% Input:
% dX      ... a vector of differentiation.
%   dX = [0 0 1] is the first derivative with respect do Z: Dz.
%   dX = [0 0 2] is the seceond derivative wrt Z: Dzz.
%   dX = [1 1 1] is Dxyz.
%
% Xcenter ... the coordinate of the (stencil) center point
% mon     ... the matrix of the exponentials obtained from the polyMat function
% dim     ... the number of dimensions.
%
% Output:
% c ... vector of the differentiated monomial terms evaluated in the center
%       point of the stencil.

import Basis.*

n = size(mon,1);
% Take the given derivative of the monomial basis.
mon_d = mon - repmat(dX, size(mon,1), 1);

% Now compute the coefficients of the differentiated monomials.
coeff = zeros(n,dim);

for l=1:dim
    dvec = repmat([0:dX(l)-1], n, 1);
    coeff(:,l) = prod(mon(:,l) - dvec, 2);
end

    
% Find the differentiatied monomial terms which due to the differentiation vanished to 0.
[idx_del,~] = find(mon_d < 0);
coeff(idx_del,:) = repmat(0, length(idx_del), dim);
mon_d(idx_del,:) = 0;

% Take the product of the coefficients dimension-wise.
coeff = prod(coeff,2)';

% Evaluate the monomials using the monomial exponential matrix and the
% coefficients
coeff_arr = repmat(coeff, size(Xcenter,1), 1);
c = coeff_arr.*polyMat_eval(Xcenter, mon_d, dim);

%
% One can also forget about all these computations if the center of the
% stencil is in (0,0,..,0). Then the coefficients are, for example in 2D:
% [0 0 1 0 0 .... 0] for the first derivative in the x-direction.
%
%coeff = zeros(n,dim);
