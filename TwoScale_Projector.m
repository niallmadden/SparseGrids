function P=TwoScale_Projector(N, sigmaN)
% TwoScale_Projector   Compute the matrix P that maps a bilinear function
%   defined on a pair of two-scale sparse grid meshes with, respectively,
%   N-by-sigmaN and (N-sigmaN)-by-sigma intervals in each direction, onto 
%   a bilinear function defined on a uniform N-by-N mesh.
% INPUT:
%   N: [1x1] integer number of intervals in the N-by-N uniform mesh
%   sigmaN: [1x1] integer proper divisor of N.
% OUTPUT:
%   P: [(N-1)^2,(2*N-N2-1)*(N2-1)] two-scale sparse grid projector matrix
%
% This function is part of 
% SparseGrids - a set of MATLAB/Octave scripts and functions
%    for implementing spare grid finite element methods.
% Authors:  Stephen Russell and Niall Madden, NUI Galway. 
% Date:     Novemberr 2016.
% Download: https://github.com/niallmadden/SparseGrids
% It is used to generate results in
%    Russell, S., and  Madden, N. An analysis and implementation of
%       sparse grid finite element methods. http://arxiv.org/abs/1511.07193
% See also Test_FEM | MultiScale_Projector

assert(mod(N,sigmaN)==0, 'Build_TwoScale_Projector(N,sigmaN): sigmaN must divide N');

x = linspace(0,1,N+1);   % 1D mesh
x2 = x(1:N/sigmaN:N+1);  % 1D sub-mesh

% The matrix p projects a mesh function on x2 onto x.
p = sparse(interp1(x2, eye(length(x2)), x));

% P1 projects from the N-by-sigmaN mesh to the N-by-N mesh
P1 = kron(p(2:end-1,2:end-1), speye(N-1));
% P2 projects from the sigmaN-by-N mesh to (N\sigmaN)-by-(N) mesh
% Identify index of points in x that are not in x2
UniqueNodes=sparse(setdiff(1:(N-1), N/sigmaN:N/sigmaN:N-N/sigmaN));
P2 = kron(sparse(UniqueNodes, 1:length(UniqueNodes),1),p(2:end-1,2:end-1));

P = [P1,P2];