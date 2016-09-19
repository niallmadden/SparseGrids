%% Test_FEM.m
% Solve the linear, two-dimensional partial differential equation
%   -(u_xx + u_yy) + u = f on (0,1)x(0,1)
% with u=0 on the boundary using
%    (1) a standard Galerkin FEM with bilinear elements,
%    (2) a two-scale sparse grid method,
%    (3) a multiscale sparse grid method,
% on a uniform mesh with N intervals in each coordinate direction.
% Authors: Niall Madden and Stephen Russell, NUI Galway. 
% Date:    September 2016.
% This code is used to generate results in
%    Russell, S., and  Madden, N. An analysis and implementation of
%       sparse grid finite element methods. http://arxiv.org/abs/1511.07193
% The code is hosted at https://github.com/niallmadden/SparseGrids
% See also FEM_System_Matrix | FEM_RHS | TwoScale_Projector | MultiScale_Projector

clear;
%% Choose which method to implement (or modify the code to use all 3)
Method = 'classical'; % Pick one of 'classical', 'two-scale', 'multiscale'

fprintf('\n\n----------------------------------------------------\n');
fprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n', ...
   'Test_FEM: solve the two-dimensional linear PDE ', ...
   '-(u_xx + u_yy) + u = f on (0,1)x(0,1) ',...
   ' with u=0 on the boundary using', ...
   '   (1) a standard Galerkin FEM with bilinear elements', ...
   '   (2) a two-scale sparse grid FEM', ...
   '   (3) a multiscale sparse grid FEM', ...
   ' on a uniform mesh with NxN grid.');

%% True solution and test Problem (true solution defined as chebfun2)
% We use a chebfun2 since it allows us to automatically compute the
% right-hand side, and a highly accurate estimate for the energy
u = chebfun2(@(x,y)4*sin(pi*x).*(y-y.^2));
f = -diffx(u,2)-diffy(u,2)+u;
% If you don't have chebfun, use this:
% u =  @(x,y)4*sin(pi*x).*(y-y.^2);
% f = @(x,y)(-4*sin(pi*x).*(pi^2*y.^2 - pi^2*y - 2) + u(x,y));
% To estimate the error, we need to compute the energy (f,u)
Energy = integral2(u.*f,[0, 1, 0, 1]);
% If you don't have chebfun, then use this:
% Energy =  5.565227840290526;

%% Discretisation data: Compute the solution for N=Ns(1), Ns(2), ...
if strcmp(Method, 'two-scale')
   Ns = (3:10).^3; % Good for two-scale
   sigma = @(n)(round(n.^(1/3)));
else
   Ns = 2.^(3:10); % Good for Galerkin and multiscale
end
%% Some matrices to store the results
Error          = nan(1,length(Ns));
Solver_Time    = nan(1,length(Ns));

%% Header for table of results
fprintf('\n---- Method: %6s ----\n', Method);
fprintf(' %4s | %9s | %6s ',...
   'N', 'Error', 'Time (s)');
fprintf('\n----------------------------\n');

%%
n=0;
for N = Ns
   n=n+1;
   %% Construct the linear system Ax=b and solve it
   A = FEM_System_Matrix(N);
   b = FEM_RHS(f, N);
   
   switch Method
      case 'classical' % Classical Galerkin method
         P = speye(length(b)); % A trivial projection for the Galerkin method
      case 'two-scale'    % the two-scale sparse grid FEM
         N2 = sigma(N);   % This should be an integer!
         P = TwoScale_Projector(N, N2);
      case 'multiscale' % the multiscale sparse grid FEM
         k=log2(N)-1;   % This should be an integer!
         P = [];
         if (mod(k,2)==1) % k is odd
            for l=0:(k-1)/2
               P = cat(2, P, MultiScale_Projector(N, N/2^l, N/2^(k-l), 2, 1));
            end
            P = cat(2, P, ...
               MultiScale_Projector(N, N/2^((k+1)/2), N/2^((k-1)/2), 1, 1));
            for l=(k+3)/2:k
               P = cat(2, P, MultiScale_Projector(N, N/2^l, N/2^(k-l), 1, 2));
            end
         else % k even
            for l=0:(k/2-1)
               P = cat(2, P, MultiScale_Projector(N, N/2^l, N/2^(k-l), 2, 1));
            end
            P = cat(2, P, MultiScale_Projector(N, N/2^(k/2), N/2^(k/2), 1, 1));
            for l=(k/2+1):k
               P = cat(2, P, MultiScale_Projector(N, N/2^l, N/2^(k-l), 1, 2));
            end
         end
      otherwise
         warning('Warning: method "%s" is not implemented', Method);
         break;
   end
   A = P'*A*P;  % The system matrix w.r.t. the new basis 
   A = (A+A')/2; % Ensure that mldivide will detect that A is symmetric

   Solver_Start=tic;
   uN = A\(P'*b);
   Solver_Time(n) = toc(Solver_Start);
   uN = P*uN;
   clear A P;
   %% Estimate the error in the energy norm
   Galerkin_Energy =  uN'*b;
   Error(n) = sqrt(Energy - Galerkin_Energy);
   clear b;
   %% Display the errors
   fprintf(' %4d | %9.3e | %7.4f \n', ...
      N, Error(n), Solver_Time(n));
end

figure(1);
loglog(Ns, Error, '-o', Ns, Ns.^(-1), '--', ...
   'LineWidth', 2, 'MarkerSize', 8);
legend('Error', 'N^{-1}');

%% % Un-comment the following to visualise results.
% [X,Y]=meshgrid(linspace(0,1,N+1));
% uN_with_BCs = zeros(N+1,N+1);
% uN_with_BCs(2:N, 2:N) = reshape(uN, (N-1), [])';
% figure(2);
% surf(X,Y,u(X,Y));
% xlabel('x'); ylabel('y');
% title('True solution');
% 
% figure(3);
% surf(X, Y, uN_with_BCs);
% xlabel('x'); ylabel('y');
% title('Finite element solution');
% 
% figure(4);
% surf(X, Y, abs(u(X,Y)-uN_with_BCs));
% xlabel('x'); ylabel('y');
% title('Error');