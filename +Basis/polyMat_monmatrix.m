function mon = polyMat_monmatrix(deg, dim)

import Basis.monomials.*

% How many monomial term does the given polynomial degree deg have?
n = mono_upto_enum(dim, deg);

% Create a matrix of the monomial exponents
mon = zeros(n,dim);
for k=2:n
    % Note: The initial monomial always has the exponent 0. This is why we
    % start at mon(k-1,:).
    mon(k,:) = mono_next_grlex(dim, mon(k-1,:));
end