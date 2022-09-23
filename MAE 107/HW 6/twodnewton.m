function [xf,yf]=twodnewton(f1,f2,x0,y0,n)
%
% Objective:
%   Find the roots of 2 functions of 2 variables)
%   
%
% input variables:
%   n - number of iterations
%   x0 - starting x
%   y0 - starting y
%   f1 - function 1
%   f2 - function 2
% output variables:
%   xf = final x of root
%   yf = final y of root
% functions called:
%   none
% Define X matrix
X = [x0;y0];
%
% Define F matrix
%
syms F(x,y);
F(x,y) = [f1;f2];
% Define Fprime matrix
%
syms Fprime(x,y);
Fprime(x,y) = [diff(f1,x) diff(f1, y) ; diff(f2,x) diff(f2,y)];
%
% For loop to apply iterations
%
for i = 1:n
    X = X - (((Fprime(X(1),X(2)))^-1)*F(X(1),X(2)));
end
xf=double(X(1));
yf=double(X(2));
%end