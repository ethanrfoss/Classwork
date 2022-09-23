function test_problem5

%
% Objective:
%   Solve a compare runge kutta solutions of different orders
% input variables:
%   none
%
% output variables:
%   none
%
% functions called:
%   fourthorderrungekutta - Solve a differential equation using fourth order runge kutta
%   secondorderrungekutta - Solve a differential equation using second order runge kutta
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
% Define fs
% 
f1= @(t,x) x/2; % x = [0 2]
f2= @(t,x) x; % x = [2 3]
f3= @(t,x) 4; % x = [3 4]
%
% Get different xk and tk vectors for different step sizes for second
% n = 10
[stk11,sxk11]=secondordererungekutta(0,2,1,f1,.5,20^1/2);
[stk12,sxk12]=secondordererungekutta(2,3,sxk11(end),f2,.5,20^1/4);
[stk13,sxk13]=secondordererungekutta(3,4,sxk12(end),f3,.5,20^1/4);
[stk1] = [stk11 stk12 stk13];
[sxk1] = [sxk11 sxk12 sxk13];
% n = 10^2
[stk21,sxk21]=secondordererungekutta(0,2,1,f1,.5,10^2/2);
[stk22,sxk22]=secondordererungekutta(2,3,sxk21(end),f2,.5,10^2/4);
[stk23,sxk23]=secondordererungekutta(3,4,sxk22(end),f3,.5,10^2/4);
[stk2] = [stk21 stk22 stk23];
[sxk2] = [sxk21 sxk22 sxk23];
% n = 10^3
[stk31,sxk31]=secondordererungekutta(0,2,1,f1,.5,10^3/2);
[stk32,sxk32]=secondordererungekutta(2,3,sxk31(end),f2,.5,10^3/4);
[stk33,sxk33]=secondordererungekutta(3,4,sxk32(end),f3,.5,10^3/4);
[stk3] = [stk31 stk32 stk33];
[sxk3] = [sxk31 sxk32 sxk33];
% n = 10^4
[stk41,sxk41]=secondordererungekutta(0,2,1,f1,.5,10^4/2);
[stk42,sxk42]=secondordererungekutta(2,3,sxk41(end),f2,.5,10^4/4);
[stk43,sxk43]=secondordererungekutta(3,4,sxk42(end),f3,.5,10^4/4);
[stk4] = [stk41 stk42 stk43];
[sxk4] = [sxk41 sxk42 sxk43];
% n = 10^5
[stk51,sxk51]=secondordererungekutta(0,2,1,f1,.5,10^5/2);
[stk52,sxk52]=secondordererungekutta(2,3,sxk51(end),f2,.5,10^5/4);
[stk53,sxk53]=secondordererungekutta(3,4,sxk52(end),f3,.5,10^5/4);
[stk5] = [stk51 stk52 stk53];
[sxk5] = [sxk51 sxk52 sxk53];
%
% Get different xk and tk vectors for different step sizes for fourth
% n = 10
[ftk11,fxk11]=fourthorderrungekutta(0,2,1,f1,20^1/2);
[ftk12,fxk12]=fourthorderrungekutta(2,3,fxk11(end),f2,20^1/4);
[ftk13,fxk13]=fourthorderrungekutta(3,4,fxk12(end),f3,20^1/4);
[ftk1] = [ftk11 ftk12 ftk13];
[fxk1] = [fxk11 fxk12 fxk13];
% n = 10^2
[ftk21,fxk21]=fourthorderrungekutta(0,2,1,f1,10^2/2);
[ftk22,fxk22]=fourthorderrungekutta(2,3,fxk21(end),f2,10^2/4);
[ftk23,fxk23]=fourthorderrungekutta(3,4,fxk22(end),f3,10^2/4);
[ftk2] = [ftk21 ftk22 ftk23];
[fxk2] = [fxk21 fxk22 fxk23];
% n = 10^3
[ftk31,fxk31]=fourthorderrungekutta(0,2,1,f1,10^3/2);
[ftk32,fxk32]=fourthorderrungekutta(2,3,fxk31(end),f2,10^3/4);
[ftk33,fxk33]=fourthorderrungekutta(3,4,fxk32(end),f3,10^3/4);
[ftk3] = [ftk31 ftk32 ftk33];
[fxk3] = [fxk31 fxk32 fxk33];
% n = 10^4
[ftk41,fxk41]=fourthorderrungekutta(0,2,1,f1,10^4/2);
[ftk42,fxk42]=fourthorderrungekutta(2,3,fxk41(end),f2,10^4/4);
[ftk43,fxk43]=fourthorderrungekutta(3,4,fxk42(end),f3,10^4/4);
[ftk4] = [ftk41 ftk42 ftk43];
[fxk4] = [fxk41 fxk42 fxk43];
% n = 10^5
[ftk51,fxk51]=fourthorderrungekutta(0,2,1,f1,10^5/2);
[ftk52,fxk52]=fourthorderrungekutta(2,3,fxk51(end),f2,10^5/4);
[ftk53,fxk53]=fourthorderrungekutta(3,4,fxk52(end),f3,10^5/4);
[ftk5] = [ftk51 ftk52 ftk53];
[fxk5] = [fxk51 fxk52 fxk53];
%
% Compute log 10 of errors for fourth and order
%
slog10e1 = log10(abs(sxk1(end)-fxk5(end)));
slog10e2 = log10(abs(sxk2(end)-fxk5(end)));
slog10e3 = log10(abs(sxk3(end)-fxk5(end)));
slog10e4 = log10(abs(sxk4(end)-fxk5(end)));
flog10e1 = log10(abs(fxk1(end)-fxk5(end)));
flog10e2 = log10(abs(fxk2(end)-fxk5(end)));
flog10e3 = log10(abs(fxk3(end)-fxk5(end)));
flog10e4 = log10(abs(fxk4(end)-fxk5(end)));
%
% Compute log of step sizes
%
log10n1 = log10(20);
log10n2 = log10(10^2);
log10n3 = log10(10^3);
log10n4 = log10(10^4);
%
% Plot solutions for different step sizes of second order
%
figure(1); hold on;
title('Solutions to diff eq for various number of steps for second order');
xlabel('t');
ylabel('x');
plot(stk1,sxk1);
plot(stk2,sxk2);
plot(stk3,sxk3);
plot(stk4,sxk4);
legend('20','10^2','10^3','10^4');
%
% Plot solutions for different step sizes of fourth order
%
figure(2); hold on;
title('Solutions to diff eq for various number of steps for fourth order');
xlabel('t');
ylabel('x');
plot(ftk1,fxk1);
plot(ftk2,fxk2);
plot(ftk3,fxk3);
plot(ftk4,fxk4);
legend('20','10^2','10^3','10^4');
%
% Plot log of step sizes vs log of errors for both
%
figure(3); hold on;
title('log10 of step sizes vs log10 of errors');
xlabel('log10 of step sizes');
ylabel('log10 of error');
plot([log10n1 log10n2 log10n3 log10n4],[slog10e1,slog10e2,slog10e3,slog10e4]);
plot([log10n1 log10n2 log10n3 log10n4],[flog10e1,flog10e2,flog10e3,flog10e4]);
legend('Second order','Fourth order');