function P=MultiScale_Projector(N, Mx, My, skip_x, skip_y)
%  MultiScale_Projector   Compute the matrix P that maps a bilinear
%   function on a uniform 2D mesh with N intervals in each coordinate
%   direction, onto a sparse grid Mx-by-My mesh, omitting each skip_x
%   and skip_y point in the Mx-mesh and My-mesh, respectively. 
% INPUT:
%   N: [1x1] integer number of x- and y-intervals in the N-by-N uniform mesh
%   Mx: [1x1] integer number of x-intervals in the (Mx)x(My) uniform mesh
%   My: [1x1] integer number of y-intervals in the (Mx)x(My) uniform mesh
%   skip_x: integer number of x-intervals to skip (should be 1 or 2)
%   skip_y: integer number of y-intervals to skip (should be 1 or 2)
% OUTPUT:
%   P:  multiscale sparse grid projector matrix 
%
% This function is part of 
% SparseGrids - a set of MATLAB/Octave scripts and functions
%    for implementing spare grid finite element methods.
% Authors:  Stephen Russell and Niall Madden, NUI Galway. 
% Date:     September 2016.
% Download: https://github.com/niallmadden/SparseGrids
% DOI:      10.5281/zenodo.154427  
% It is used to generate results in
%    Russell, S., and  Madden, N. An analysis and implementation of
%       sparse grid finite element methods. http://arxiv.org/abs/1511.07193

% Authors: Niall Madden and Stephen Russell, NUI Galway. 
% Date:    September 2016.
% This code is used to generate results in
%   Russell, S., and  Madden, N. An analysis and implementation of
%     sparse grid finite element methods. http://arxiv.org/abs/1511.07193
% The code is hosted at https://github.com/niallmadden/SparseGrids
% See also Test_FEM | TwoScale_Projector

assert(mod(N,Mx)==0, 'MultiScale_Projector: Mx must divide N');
assert(mod(N,My)==0, 'MultiScale_Projector: My must divide N');
assert(skip_x==1 || skip_x==2, 'MultiScale_Projector: skip_x must be 1 or 2');
assert(skip_y==1 || skip_y==2, 'MultiScale_Projector: skip_y must be 1 or 2');

x = linspace(0,1,N+1); % 1D mesh

Sx = x(1:N/Mx:N+1);     % 1D sparse mesh
Sy = x(1:N/My:N+1);     % 1D sparse mesh

% The matrices px and py project from the Sx- and Sy-meshes to the
% x- and y-meshes respectively
px = sparse( interp1(Sx, eye(length(Sx)), x));
py = sparse( interp1(Sy, eye(length(Sy)), x));

P = kron(py(2:N, 2:skip_y:My),px(2:N, 2:skip_x:Mx));