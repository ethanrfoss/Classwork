function test_fourthorderrungekutta

%
% Objective:
%   Solve a compare runge kutta solutions of different step sizes
% input variables:
%   none
%
% output variables:
%   none
%
% functions called:
%   fourthorderrungekutta - Solve a differential equation using fourth order runge kutta
%
%
%
% Erase all old variables in workspace.
%
clear;
%
% Close all the exiting plot windows.
%
close all;
%
% Clear the desktop.
%
clc;
% 
% Define f
% 
f= @(t,x) t + (sin(t)*(x^2 + 1)^(1/2))/2;
%
% Get different xk and tk vectors for different step sizes
%
[tk1,xk1]=fourthorderrungekutta(0,5,1,f,10^1);
[tk2,xk2]=fourthorderrungekutta(0,5,1,f,10^2);
[tk3,xk3]=fourthorderrungekutta(0,5,1,f,10^3);
[tk4,xk4]=fourthorderrungekutta(0,5,1,f,10^4);
[tktrue,xktrue]=fourthorderrungekutta(0,5,1,f,10^5);
%
% Compute log 10 of errors
%
log10e1 = log10(abs(xk1(end)-xktrue(end)));
log10e2 = log10(abs(xk2(end)-xktrue(end)));
log10e3 = log10(abs(xk3(end)-xktrue(end)));
log10e4 = log10(abs(xk4(end)-xktrue(end)));
%
% Compute log of step sizes
%
log10n1 = log10(10);
log10n2 = log10(10^2);
log10n3 = log10(10^3);
log10n4 = log10(10^4);
%
% Plot solutions for different step sizes
%
figure(1); hold on;
title('Solutions to diff eq for various number of steps');
xlabel('t');
ylabel('x');
plot(tk1,xk1);
plot(tk2,xk2);
plot(tk3,xk3);
plot(tk4,xk4);
legend('10','10^2','10^3','10^4');
%
% Plot log of step sizes vs log of errors
%
figure(2); hold on;
title('log10 of step sizes vs log10 of errors');
xlabel('log10 of step sizes');
ylabel('log10 of error');
plot([log10n1 log10n2 log10n3 log10n4],[log10e1,log10e2,log10e3,log10e4]);