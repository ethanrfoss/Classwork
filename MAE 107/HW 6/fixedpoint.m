function [xk,flag,diverged]=fixedpoint(g,x0,n, minerror)
%
% Objective:
%   Find the roots of a function using fixed point method
%   
%
% input variables:
%   n - max number of iterations
%   x0 - starting x
%   g - function to find root of, adjusted for fixed point
%   minerror - minimum error necessary to break loop
% output variables:
%   xk - iterates
%   flag - string displaying whether max steps were reached or
%   not
% functions called:
%   none
% Define length of matrix and intitial condition
xk = zeros(n,1);
xk(1) = x0;
%
% For loop to apply iterations
%
flag =false;
diverged = false;
for i = 1:(n-1)
    xk(i+1) = g(xk(i));
    if abs(xk(i+1)-xk(i))<=minerror
        flag = true;
        xk((i+2):end) = [];
        break;
    end
    if xk(i+1) == inf
        diverged = true;
        xk((i+2):end) = [];
        break;
    end
end
%end