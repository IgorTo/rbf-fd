function wghts = generateWeights(X_stencil, Y_center, idx_Y_X, idx, Minv, x_center, basis, epOrP, polydeg, n, dim, polys)
import RBFEngine.Basis.*


    % Take the interpolation matrix from x in X which is closest to the current y in Y.
    idx_inv = idx_Y_X(k);    
    
    % The Y-center.
    idx_c = k;
    
    %
    % Shift and scale the stencil nodes to a circle-unit.
    %
    Xscaled = zeros(length(idx),dim);
    Yscaled = zeros(1,dim);
    
    % Shift the stencil nodes so that the center node is the origin.
    
    Xscaled = X_stencil - repmat(X_stencil(1,:), n, 1); %
    
    % Do the same for the Y-center.
    Yscaled = Y_center - repmat(X_stencil(1,:), n, 1); %

    % Scale the points to [0,1]x[0,1]x...x[0,1]
    scale = 1./max(abs(X_local));
    scale_stencil = repmat(scale, n, 1);

    Xscaled = Xscaled .* scale_stencil;
    Yscaled = Yscaled .*scale_stencil(1,:);
    
    %
    % Stencil distance matrices. Center point is in the origin (0).
    %
    dx_stencil = Yscaled(1) - Xscaled(:,1);
    dy_stencil = Yscaled(2) - Xscaled(:,2);

    %
    % Distance matrices.
    %
    dX = cell(dim,1); % first cell x-x', second cell y-y',...
    for j=1:dim
        dX{j} = Yscaled(j) - Xscaled(:,j);
    end
    
    
%     % Put the differences to eps to avoid divisions by 0.
%     if (dx_stencil(1) == 0)
%         dx_stencil(1) = eps;
%     end
%     
%     if(dy_stencil(1) == 0)
%         dy_stencil(1) = eps;
%     end
    
switch dim
    case 1
        r_stencil = r(dX{1});
    case 2
        r_stencil = r(dX{1}, dX{2});
    case 3
        r_stencil = r(dX{1}, dX{2}, dX{3});
    case 4
        r_stencil = r(dX{1}, dX{2}, dX{3}, dX{4});
end
    
    %
    % 0th derivatives
    %
    b = bf(r_stencil,epOrP); % Here the distance is actually not scaled.        

    %
    % 1st derivatives
    %    
    bx = Dbf(r_stencil, epOrP, dX{1});
    by = Dbf(r_stencil, epOrP, dX{2});

    %
    % Second derivatives
    %
    bxx = DDbf(r_stencil, epOrP, dX{1});
    byy = DDbf(r_stencil, epOrP, dX{2});
    bxy = DDbf_ij(r_stencil, epOrP, dX{1}, dX{2});
    
    %
    % Compute the stencils at once.
    %
    stenc =  Minv * [bx by bxx byy bxy b; polys.cx(k,:)' polys.cy(k,:)' polys.cxx(k,:)' polys.cyy(k,:)' polys.cxy(k,:)' polys.c(k,:)'];
    
    wghts.Dx_loc(k, :) = scale(1).*stenc(1:idx_st(1)-n_mon);
    wghts.Dy_loc(k, :) = scale(2).*stenc(idx_st(1)+1:idx_st(2)-n_mon);
    wghts.Dxx_loc(k, :) = scale(1)^2.*stenc(idx_st(2)+1:idx_st(3)-n_mon);
    wghts.Dyy_loc(k, :) = scale(2)^2.*stenc(idx_st(3)+1:idx_st(4)-n_mon);
    wghts.Dxy_loc(k, :) = scale(1)*scale(2).*stenc(idx_st(4)+1:idx_st(5)-n_mon);
    wghts.Interp_loc(k,:) = stenc(idx_st(5)+1:idx_st(6)-n_mon);
end