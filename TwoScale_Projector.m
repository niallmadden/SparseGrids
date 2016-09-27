function P=TwoScale_Projector(N, N2)
% TwoScale_Projector   Compute the matrix P that maps a bilinear function
%   defined on a pair of two-scale sparse grid meshes with, respectively,
%   N-by-N2 and (N-N2)-by-N2 intervals in each direction, onto a bilinear
%   function defined on a uniform N-by-N mesh.
% INPUT:
%   N: [1x1] integer number of intervals in the N-by-N uniform mesh
%   N2: [1x1] integer proper divisor of N.
% OUTPUT:
%   P: [(N-1)^2,(2*N-N2-1)*(N2-1)] two-scale sparse grid projector matrix
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
% See also Test_FEM | MultiScale_Projector

assert(mod(N,N2)==0, 'Build_TwoScale_Projector(N,N2): N2 must divide N');

x = linspace(0,1,N+1);   % 1D mesh
x2 = x(1:N/N2:N+1);      % 1D sub-mesh

% The matrix p projects from the N2-mesh to the N-mesh
p = sparse(interp1(x2, eye(length(x2)), x));

% P1 projects from the N-by-N2 mesh to N-by-N mesh
P1 = kron(p(2:end-1,2:end-1), speye(N-1));
% P2 projects from the N2-by-N mesh to (N/N2)-by-(N) mesh
% Identify index of points in x that are not in x2
UniqueNodes=sparse(setdiff(1:(N-1), N/N2:N/N2:N-N/N2));
P2 = kron(sparse(UniqueNodes, 1:length(UniqueNodes),1),p(2:end-1,2:end-1));

P = [P1,P2];