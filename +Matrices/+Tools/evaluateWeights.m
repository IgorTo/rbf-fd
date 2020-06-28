function b_all = evaluateDerivatives(X, Y, basis, epOrP)
% This function evaluates the derivatives in the order:
% [bx by bxx byy bxy b]

dim = size(X,2);
%
% Parse the basis functions.
%
Dbf = basis.Dbf;
DDbf = basis.DDbf;
DDbf_ij = basis.DDbf_ij;
bf = basis.bf;
r = basis.r;

[r_stencil, dX] = RBFEngine.LeastSquares.generateDistanceMatrix(X, Y, r);

%
% 0th derivatives
%
b = bf(r_stencil,epOrP); % Here the distance is actually not scaled.

%
% 1st derivatives
%
b1 = cell(dim,1);
for k=1:dim
    b1{k} = Dbf(r_stencil, epOrP, dX{k});
end

%
% Second derivatives
%
n_combs = nchoosek(dim, 2);
b2 = cell(dim + n_combs,1);

for k=1:dim
    b2{k} = DDbf(r_stencil, epOrP, dX{k});
end

% Mixed 2nd derivatives.
combs = nchoosek(1:dim, 2);
for k=1:size(combs,1)
    idx1 = combs(k,1);
    idx2 = combs(k,2);
    
    b2{k} = DDbf_ij(r_stencil, epOrP, dX{idx1}, dX{idx2});
end

% Join all of the derivatives together.
b_all = [cell2mat(b1); cell2mat(b2); b]
