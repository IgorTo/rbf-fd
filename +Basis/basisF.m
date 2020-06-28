function basis = basisF(bf, dim)



%
% All the basis functions should be written in terms of the distance
% function, the shape parameter or the PHS degree, and the independent
% variable(s), depending on the order of the derivative.
%

%
% Now choose the basis function.
%
switch bf
    %
    % PHS (Polyharmonic spline) basis
    %
    case 'PHS'
        basis.bf = @(r, p) r.^p;
        basis.Dbf = @(r,p,x) p.*x.*r.^(p-2);
        basis.DDbf = @(r,p,x) p.*r.^(p-2) + p*(p-2).*r.^(p-4).*x.^2;
        basis.DDbf_ij = @(r,p,x,y) p*(p-2)*r.^(p-4).*x.*y;
        
    case 'Gaussian'
        basis.bf = @(r, p) exp(-(p^2.*r.^2));
        basis.Dbf = @(r, p, x) exp(-(p^2.*r.^2)) .* (-2 * p^2 .*x); % derivative in x-direction
        basis.DDbf = @(r, p, x) exp(-(p^2.*r.^2)) .* (-2.*p^2 + 4.* p^4 .* x.^2); % derivative in x-direction
        basis.DDbf_ij = @(r, p, x,y) exp(-(p^2.*r.^2)) * 4.* p^4 .* x .* y;

    case 'MQ'
        basis.bf = @(r, p) (1 + p^2.*r.^2).^(1/2);
        basis.Dbf = @(r, p, x) p^2./((1 + p^2.*r.^2).^(1/2)) .* x;
        basis.DDbf = @(r, p, x) p^2./((1 + p^2.*r.^2).^(1/2)) - p^4./((1 + p^2.*r.^2).^(3/2)) .* x.^2;
        basis.DDbf_ij = @(r,p,x,y) p^4./((1 + p^2.*r.^2).^(3/2)) .* x .* y;

end


%
% Choose the dimension. Only the distance function changes.
%
switch dim
    case '1d'
        basis.r = @(x) sqrt(x.^2);
    case '2d'
        basis.r = @(x,y) sqrt(x.^2 + y.^2);
    case '3d'
        basis.r = @(x,y,z) sqrt(x.^2 + y.^2 + z.^2);
end