# RBF-FD matrices
The core of this slim repository is a Matlab script that constructs RBF-FD differentiation and evaluation matrices to use in least-squares and collocation PDE methods.

## Usage


```matlab
[E, Dx, Dy, Dxx, Dyy, Dxy, stencils] = Matrices.generate_2d(X, Y, p, n, polydeg, bf, indeces, []),
```
Input:
X ... Nx1 list of interpolation points.
Y ... Mx1 list of evaluation points. (If Y=X, then the matrices will be square)
p ... exponent on PHS
n ... stencil size.
polydeg ... polynomial degree.
basis ... the basis functions.
indeces ... a struct with indeces of stencil neighbors (idx_X) and indeces of closest stencils to every evaluation points (idx_Y_X)
stencils ... a struct of already computed stencils (can be used for interpolation purposes). Input [] when no precomputed stencils are available.

Output:
E ... MxN evaluation matrix.
Dx ... MxN differentiation matrix d/dx
Dy ... MxN differentiation matrix d/dy
Dxx ... MxN differentiation matrix d^2/dx^2
Dxx ... MxN differentiation matrix d^2/dy^2
Dxy ... MxN differentiation matrix d^2/(dx*dy)


