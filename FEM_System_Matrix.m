function A = FEM_System_Matrix(N)
% FEM_System_Matrix   Compute the system matrix for the standard
%   Galerkin FEM with bilinear elements applied to
%   (u_x,v_x) + (u_y,v_y) + (u,v) = (f,v)  on  the unit square, 
%   on a uniform tensor product mesh with N intervals in each coordinate
%   direction.
% INPUT:
%   N: [1x1] integer number of intervals for the NxN uniform mesh.
% OUTPUT:
%   A: [(N-1)^2,(N-1)^2] sparse system matrix.
% 
% This function is part of 
% SparseGrids - a set of MATLAB/Octave scripts and functions
%    for implementing spare grid finite element methods.
% Authors:  Stephen Russell and Niall Madden, NUI Galway. 
% Date:     September 2016.
% Download: https://github.com/niallmadden/SparseGrids
% It is used to generate results in
%    Russell, S., and  Madden, N. An analysis and implementation of
%       sparse grid finite element methods. http://arxiv.org/abs/1511.07193
% See also Test_FEM | FEM_RHS 

%%  1D stencils, and associated linear systems
s0 = 1/(6*N)*[1,4,1];
s2 = N*[-1,2,-1];
a0 = StencilToTridiagonal(s0, N-1); % 1 dimensional mass matrix
a2 = StencilToTridiagonal(s2, N-1); % 1 dimensional stiffness matrix

%% The 2D system matrix
% One can express the system matrix as
% A  = kron(a0,a2) + kron(a2,a0) + kron(a0,a0);
% where the three products come from, respectively, the terms
%  (u_x,v_x), (u_y,v_y) and from (u,v) in the bilinear form.
% But it is more efficient to compute this as:
A = kron(a0, a2) + kron(a2+a0,a0);
end


function A=StencilToTridiagonal(s, n)
% Form an (n)x(n) tridiagonal matrix from a stencil
% INPUTS:
%  s:  [1x3] stencil [s(1),s(2),s(3)]=[A(i,i-1, A(i,i), A(i,i+1)]
%  n:  [1x1] number of rows and columns in the matrix.
% OUTPUT:
%  A:  [(n)x(n)] Tridiagonal matrix
A =   sparse(2:n,   1:n-1, s(1), n, n) ...
    + sparse(1:n,   1:n,   s(2)) ... 
    + sparse(1:n-1, 2:n,   s(3), n, n);
end