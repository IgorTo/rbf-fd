function [X, scale_stencil] = scaleStencil(X, x_center, varargin)

n = size(X,1);
% Shift the template so that x_center is in 0. Scale it to
% [-1 1] x [-1 1] x ... x [-1, 1]
X = X - repmat(x_center, n, 1); %

if(length(varargin)==1)
    scale = varargin{1};
else
    scale = 1./max(abs(X));
end

scale_stencil = repmat(scale, n, 1);
% Scale the stencil points.
X = X .* scale_stencil;