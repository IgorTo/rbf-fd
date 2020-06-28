function [M, M_inv] = generateInterpMatrix(X_local, basis, monMatrix, epOrP, dim)
%import([basis '.*']);
import Basis.*
%
% This function throws out an interpolation matrix.
% Input:
% X_local ... the pointset in which the matrix is to be made.
% x_center ... the stencil center point.
% basis ... an object containing basis.bf and basis.r.
% epOrP .. the shape parameter or the PHS degree
% polydeg ... the polynomial degree to be included
% n ... the stencil size
% dim ... the number of dimensions


%bf = basis.bf;
%r = basis.r;



%
% Distance matrices.
%
r_matrix = Matrices.Tools.generateDistanceMatrix(X_local, X_local);

%
% Interp. matrix
%
A_loc = basis.bf(r_matrix, epOrP);

%
% Polynomial matrix
%
P = polyMat_eval(X_local, monMatrix, dim);

%
% The RBF-system matrix with polynomial augmentation.
%
zrs = zeros(size(P,2), size(P,2));

M = [A_loc, P; ...
    P', zrs];

M_inv = inv(M);
end