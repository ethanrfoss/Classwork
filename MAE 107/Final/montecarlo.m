function [tk,xk]= montecarlo(t0,tf,x0,f,sigma,n)
%
% Objective:
%   Solve a differential equation using monte carlo method
%   
%
% input variables:
%   t0 - t initial
%   tf - t final
%   x0 - x initial
%   f(x,t) - function of derivative of x with respect to t
%   n - number of steps
%
% output variables:
%   xk = vector solution of x variables
%   tk = vector of t values
% functions called:
%   none
%
% Establish step size and vector sizes
%
h = (tf-t0)/n;
tk = t0:h:tf;
xk = length(tk);
xk(1) = x0;
% for loop to solve runge kutta
%
for i = 1:length(tk)-1
    xk(i+1) = xk(i) + h*f(tk(i),xk(i))+sigma*sqrt(h)*randn;
end
%
% end
%