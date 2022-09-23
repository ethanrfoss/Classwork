function [xn,un]=diffeqsolver(n, u0,uf)
%
% Objective:
%   Solve a diff eq using a tridiagonal matrix
%   -2*n^2*u(x(k-1))+(4*n^2+1)*u(x(k))-2*n^2*u(x(k+1)) = 61(1-(x(k)-.5)^2)
%   
%
% input variables:
%   n - number of steps
%   u0 - u initial
%   uf - u final
% output variables:
%   xn = x vector
%   un - solution vector
% functions called:
%   tridiagonal
%
%
% Create x vector for x values
%
xn = 0:1/n:1;
% Establish u,d,l to be solved using tridiagonal
%
u = -2*n^2*ones(n-2,1);
un = zeros(n+1,1);
d = (4*n^2+1)*ones(n-1,1);
l = -2*n^2*ones(n-2,1);
b = ones(n-1,1);
un(1) = u0; 
un(end) = uf;
%
% For loop to determine b matrix
%
b(1) = 61*(1-(xn(2)-.5)^2) + 2*n^2*un(1);
for i=2:n-2
    b(i) = 61*(1-(xn(i+1)-.5)^2);
end
b(end) = 61*(1-(xn(end-1)-.5)^2) + 2*n^2*un(end);
%
% Input into tridiagonal to get solution for u
%
un(2:n) = tridiagonal(u,d,l,b);