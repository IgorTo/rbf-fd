# RBF-FD matrices
The core of this slim repository is a Matlab script that constructs RBF-FD (phs+poly) differentiation and evaluation matrices which can be used in least-squares and collocation methods for solving PDEs in two dimensions.

## Tags
- RBF-FD
- PHS+poly
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
p ... exponent on PHS
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
