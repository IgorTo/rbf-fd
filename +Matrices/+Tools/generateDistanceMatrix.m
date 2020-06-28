function [r_matrix, dX] = generateDistanceMatrix(X_local, Y_local)

%
% Distance matrices.
%
dim = size(X_local,2);
n_X = size(X_local,1);
n_Y = size(Y_local,1);

dX = cell(dim,1); % first cell x-x', second cell y-y',...
for j=1:dim    
    dx_loc = repmat(X_local(:,j), 1, n_Y);
    dy_loc = repmat(Y_local(:,j), 1, n_X);
    
    if(n_Y == 1)
        dX{j} = dy_loc - dx_loc;
    else
        dX{j} = dx_loc - dy_loc';
    end
end

switch dim
    case 1
        r_matrix = sqrt(dX{1}.^2);
    case 2
        r_matrix = sqrt(dX{1}.^2 + dX{2}.^2);
    case 3
        r_matrix = sqrt(dX{1}.^2 + dX{2}.^2 + dX{3}.^2);
end