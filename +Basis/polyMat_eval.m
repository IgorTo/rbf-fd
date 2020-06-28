function P = polyMat_eval (X, mon, dim)

N = size(X,1);
n = size(mon,1);
P = ones(N,n);

XX = cell(dim,1);
mon_Tr = mon';

for d=1:dim
    XX{d} = repmat(X(:,d), 1, n);
    mmon = repmat(mon_Tr(d,:), N, 1);
    P = P .* XX{d}.^mmon;
end