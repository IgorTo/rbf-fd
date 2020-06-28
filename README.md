# RBF-FD matrices
The core of this slim repository is a Matlab script that constructs RBF-FD polynomially augmented differentiation and evaluation matrices which can be used in least-squares and collocation methods for solving PDEs in two dimensions.

## Tags
- RBF-FD
- PHS+poly, GA+poly, MQ+poly
- Collocation
- Least-squares
- Oversampled operators
- 2D

## Usage

```matlab
[E, Dx, Dy, Dxx, Dyy, Dxy, stencils] = Matrices.generate_2d(X, Y, p, n, polydeg, bf, indeces, []),
```
### Inputs
```
X ... Nx2 list of interpolation points.
Y ... Mx2 list of evaluation points. (If Y=X, then the matrices will be square)
p ... exponent on the polyharmonic spline (PHS) or the shape parameter of a gaussian or a multiquadric.
n ... stencil size.
polydeg ... polynomial degree.
basis ... the basis functions.
indeces ... a struct with indeces of stencil neighbors (idx_X) and indeces of closest stencils to every evaluation points (idx_Y_X)
stencils ... a struct of already computed stencils (can be used for interpolation purposes). Input [] when no precomputed stencils are available.
```
### Outputs
```
E ... MxN evaluation matrix.
Dx ... MxN differentiation matrix d/dx
Dy ... MxN differentiation matrix d/dy
Dxx ... MxN differentiation matrix d^2/dx^2
Dxx ... MxN differentiation matrix d^2/dy^2
Dxy ... MxN differentiation matrix d^2/(dx*dy)
```
## Example
An example of using rectangular matrices is provided in file Example_Poisson_fitted.m.

## Other packages used
The only other package that this library uses is MONOMIAL (https://people.sc.fsu.edu/~jburkardt/m_src/monomial/monomial.html) maintained under LGPL license.

# Papers
The papers which use this script are:
- I. Tominec, E. Larsson, A. Heryudono: A least squares radial basis function finite difference method with improved stability properties
- I. Tominec, E. Breznik: Unfitted radial basis function finite difference method in a least-squares setting
