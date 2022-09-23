function test_montecarlo

%
% Objective:
%   Create histogram of monte carlo method
% input variables:
%   none
%
% output variables:
%   none
%
% functions called:
%   montecarlo - Solve a differential equation using monte carlo
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
f = @(t,x) 1-sin(pi*x/4)
%
% Define size of x vectors and Use for loop to run 500 iterations twice
%
x5001 = zeros(500,1);
x5002 = zeros(500,1);
for i = 1:500
    [tk,xk] = montecarlo(0,10,1,f,.1,1000);
    x5001(i) = xk(end);
end
for i = 1:500
    [tk,xk] = montecarlo(0,10,1,f,.1,1000);
    x5002(i) = xk(end);
end
%
% Define size of x vectors and Use for loop to run 500 iterations twice
%
x20001 = zeros(1,2000);
x20002 = zeros(1,2000);
for i = 1:2000
    [tk,xk] = montecarlo(0,10,1,f,.1,1000);
    x20001(i) = xk(end);
end
for i = 1:2000
    [tk,xk] = montecarlo(0,10,1,f,.1,1000);
    x20002(i) = xk(end);
end
%
% Plot histograms
%
figure(1); hold on;
tiledlayout(2,2);
nexttile; histogram(double(x5001),20); title('Nr = 500'); xlabel('xend'); ylabel('number');
nexttile; histogram(double(x5002),20); title('Nr = 500'); xlabel('xend'); ylabel('number');
nexttile; histogram(double(x20001),20); title('Nr = 2000'); xlabel('xend'); ylabel('number');
nexttile; histogram(double(x20002),20); title('Nr = 2000'); xlabel('xend'); ylabel('number');
%
% end
%