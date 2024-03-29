# Rectangular and square RBF-FD matrices in MATLAB
The core of this slim repository is a Matlab script that constructs RBF-FD differentiation and evaluation matrices. Those can be used in the least-squares and the collocation methods for solving PDEs in two dimensions.

## Tags
- RBF-FD
- PHS+poly, GA+poly, MQ+poly
- Collocation
- Least-squares
- Oversampled operators
- 2D

## Usage

```matlab
% Parameters.
load X-pts and Y-pts
bf = Basis.basisF('PHS', '2d'); % Choose basis functions.
p = 3; % PHS power (r^p).
polydeg = 3; % Augmented polynomial degree.
n = 2*nchoosek(polydeg+2,2); % Stencil size.

% Get matrices.
[E, Dx, Dy, Dxx, Dyy, Dxy, stencils] = Matrices.generate_2d(X, Y, p, n, polydeg, bf, stencils[optional]),
```
### Inputs
```
X ... Nx2 list of interpolation points.
Y ... Mx2 list of evaluation points. (If Y=X, then the matrices will be square)
p ... exponent on the polyharmonic spline (PHS) or the shape parameter of a gaussian or a multiquadric.
n ... stencil size.
polydeg ... polynomial degree.
basis ... the basis functions: 'PHS', 'GA', 'MQ'.
indeces ... a struct with indeces of stencil neighbors (indeces.idx_X) and indeces of closest stencils to every evaluation points (indeces.idx_Y_X)
stencils (optional) ... a struct of already computed stencils (can be used for interpolation purposes). Input [] when no precomputed stencils are available.
```
### Outputs
```
E ... MxN evaluation matrix.
Dx ... MxN differentiation matrix d/dx
Dy ... MxN differentiation matrix d/dy
Dxx ... MxN differentiation matrix d^2/dx^2
Dyy ... MxN differentiation matrix d^2/dy^2
Dxy ... MxN differentiation matrix d^2/(dx*dy)
```
## Examples
1. Using the fitted RBF-FD method to solve a Poisson equation. [Example_Poisson_fitted.m](https://github.com/IgorTo/rbf-fd/blob/master/Example_Poisson_fitted.m).
2. Using the unfitted RBF-FD method to solve a Poisson equation. [Example_Poisson_unfitted.m](https://github.com/IgorTo/rbf-fd/blob/master/Example_Poisson_unfitted.m).

## Other packages used
The only other package that this library uses is MONOMIAL (https://people.sc.fsu.edu/~jburkardt/m_src/monomial/monomial.html) maintained under LGPL license.

# Papers
The papers which use this script are:
- Stability estimates for radial basis function methods applied to time-dependent hyperbolic PDEs. Preprint, 2021.
- Residual viscosity stabilized RBF-FD methods for solving nonlinear conservation laws. Preprint, 2021.
- An unfitted radial basis function generated finite difference method applied to thoracic diaphragm simulations. Preprint, 2021.
- A least squares radial basis function finite difference method with improved stability properties. SIAM Journal on Scientific Computing, 2021.
- An unfitted RBF-FD method in a least-squares setting for elliptic PDEs on complex geometries. Journal of Computational Physics, 2021.
