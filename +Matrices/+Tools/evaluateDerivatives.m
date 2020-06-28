function [b, bx, by, bxx, byy, bxy] = evaluateDerivatives(X, Y, basis, epOrP)
% This function evaluates the derivatives in the order:
% [bx by bxx byy bxy b]

dx_stencil = Y(:,1) - X(:,1);
dy_stencil = Y(:,2) - X(:,2);


r_stencil = basis.r(dx_stencil, dy_stencil); % here the distance matrix is actually not scaled.



b = basis.bf(r_stencil,epOrP); % Here the distance is actually not scaled.

bx = basis.Dbf(r_stencil, epOrP, dx_stencil);
by = basis.Dbf(r_stencil, epOrP, dy_stencil);

bxx = basis.DDbf(r_stencil, epOrP, dx_stencil);
byy = basis.DDbf(r_stencil, epOrP, dy_stencil);
bxy = basis.DDbf_ij(r_stencil, epOrP, dx_stencil, dy_stencil);

% dim = size(X,2);
% %
% % Parse the basis functions.
% %
% Dbf = basis.Dbf;
% DDbf = basis.DDbf;
% DDbf_ij = basis.DDbf_ij;
% bf = basis.bf;
% r = basis.r;
% 
% [r_stencil, dX] = RBFEngine.LeastSquares.generateDistanceMatrix(X, Y, r);
% 
% b_all = [];
% 
% %
% % 0th derivatives
% %
% b = bf(r_stencil,epOrP); % Here the distance is actually not scaled.
% b_all = [b_all, b];
% 
% 
% %
% % 1st derivatives
% %
% %b1 = cell(1, dim);
% for k=1:dim
%     %b1{k} = Dbf(r_stencil, epOrP, dX{k});
%     b_all = [b_all, Dbf(r_stencil, epOrP, dX{k})];
% end
% 
% %
% % Second derivatives
% %
% %n_combs = nchoosek(dim, 2);
% %b2 = cell(1, dim + n_combs);
% 
% for k=1:dim
%     %b2{k} = DDbf(r_stencil, epOrP, dX{k});
%     b_all = [b_all, DDbf(r_stencil, epOrP, dX{k})];
%     
% end
% 
% 
% % Mixed 2nd derivatives.
% combs = nchoosek(1:dim, 2);
% for k=1:size(combs,1)
%     idx1 = combs(k,1);
%     idx2 = combs(k,2);
%     
%     %b2{idx1+idx2} = DDbf_ij(r_stencil, epOrP, dX{idx1}, dX{idx2});
%     b_all = [b_all DDbf_ij(r_stencil, epOrP, dX{idx1}, dX{idx2})];
% end
% % Join all of the derivatives together.
% %b_all = [cell2mat(b1), cell2mat(b2), b];
