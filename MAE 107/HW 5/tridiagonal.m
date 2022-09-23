function [x]=tridiagonal(u,d,l,b)
%
% Objective:
%   Solve a tridiagonal system
%   
%
% input variables:
%   u - top diagonal
%   d - middle diagonal
%   l - bottom diagonal
%   b - b vector
%
% output variables:
%   x = solution of the augmented matrix
% functions called:
%   none
%
% Establish length of x and final u,d,b
%
x = zeros(length(d),1);
un = zeros(length(u),1);
dn = zeros(length(d),1);
bn = zeros(length(d),1);
%
% Establish initial values for un,dn,bn
%
un = u;
dn(1) = d(1);
bn(1) = b(1);
% for loop to simplify tridiagonal matrix
%
for i = 2:length(b)
    dn(i) = d(i) -(l(i-1)*un(i-1))/dn(i-1);
    bn(i) = b(i) -(l(i-1)*bn(i-1))/dn(i-1);
end
%
% For loop for back substitution
%
x(end) = bn(end)/dn(end);
for i = (length(x)-1):-1:1
   x(i) = (bn(i) - un(i)*x(i+1))/dn(i);
end
%
% 
%