function [tk,xk]=fourthorderrungekutta(t0,tf,x0,f,n)
%
% Objective:
%   Solve a differential equation using fourth order runge kutta
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
halfh = h/2;
for i = 1:length(tk)-1
    K1 = h*f(tk(i),xk(i));
    K2 = h*f(tk(i)+halfh,xk(i)+K1/2);
    K3 = h*f(tk(i)+halfh,xk(i)+K2/2);
    K4 = h*f(tk(i)+h,xk(i)+K3);
    xk(i+1) = xk(i) + K1/6+K2/3+K3/3+K4/6;
end
%
% end
%