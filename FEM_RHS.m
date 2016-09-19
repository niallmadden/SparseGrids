function b = FEM_RHS(f, N)
% FEM_RHS    Compute the right-hand side of a linear system for the
%   Galerkin FEM with bilinear elements applied to
%   (u_x,v_x) + (u_y,v_y) + (u,v) = (f,v)  on  the unit square, 
%   on a uniform tensor product mesh with N intervals in each coordinate
%   direction.
% INPUTS:
%   f: [function handle @(x,y)] right-hand side of the differential equation
%   N: [1x1] integer number of intervals for the NxN uniform mesh
% OUTPUT:
%   b: [(N-1)^2,1] vector for the right-hand side. This is used when 
%               solving A*uN=b where A is computed by FEM_System_Matrix.m
% Authors: Niall Madden and Stephen Russell, NUI Galway. 
% Date:    September 2016.
% This code is used to generate results in
%    Russell, S., and  Madden, N. An analysis and implementation of
%       sparse grid finite element methods. http://arxiv.org/abs/1511.07193
% The code is hosted at https://github.com/niallmadden/SparseGrids
% See also Test_FEM | FEM_System_Matrix

h = 1/N;
x = linspace(0,1,N+1);
[X,Y] = meshgrid(x);

%%  Define the quadrature rule: we use 2-point Gaussian Quadrature
qN = 2; % Number of points in each coordinate direction
qW = [1 1;1 1]/4;  % Weights as a qN-times-qN matrix
qx = [1-1/sqrt(3), 1+1/sqrt(3)]/2; % x-values for function evaluation
qy = [1-1/sqrt(3), 1+1/sqrt(3)]/2; % y-values for function evaluation

%% The four bilinear functions that contribute to an individual basis
%   function, evaluated at the quadrature points.
phi1 = (1-qx)'*(1-qy);   % phi1(0,0)=1
phi2 =     qx'*(1-qy);   % phi2(1,0)=1
phi3 =     qx'*qy;       % phi3(1,1)=1
phi4 = (1-qx)'*qy;       % phi4(0,1)=1

%% Assemble the right-hand side vector 
b = zeros(N-1);
for i=1:qN
   for j=1:qN
      b = b + h^2*qW(i,j)*(...
         f(X(2:N,  2:N)+h*qx(i),   Y(2:N,  2:N)+h*qy(j))*phi1(i,j) ...
         + f(X(2:N,  1:N-1)+h*qx(i), Y(2:N,  1:N-1)+h*qy(j))*phi2(i,j) ...
         + f(X(1:N-1,1:N-1)+h*qx(i), Y(1:N-1,1:N-1)+h*qy(j))*phi3(i,j) ...
         + f(X(1:N-1,2:N)+h*qx(i),   Y(1:N-1,2:N)+h*qy(j))*phi4(i,j));
   end
end
b = reshape(b', [], 1);