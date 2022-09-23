function [tk,xk]=secondorderrungekutta(t0,tf,x0,f,gamma,n)
%
% Objective:
%   Solve a differential equation using second order runge kutta
%   
%
% input variables:
%   t0 - t initial
%   tf - t final
%   x0 - x initial
%   f(x,t) - function of derivative of x with respect to t
%   gamma - gamma for runge kutta
%   n - number of steps
%
% output variables:
%   xk = vector solution of x variables
%   tk = vector of t values
% functions called:
%   none
%
% Establish beta, alpha and upsilon
%
beta = 1 - gamma;
alpha = .5/gamma;
upsilon = .5/gamma;
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
     K1 = h*f(tk(i),xk(i));
     K2 = h*f(tk(i)+upsilon*h,xk(i)+alpha*K1);
     xk(i+1) = xk(i) + beta*K1+gamma*K2;
 end
%
%
% end