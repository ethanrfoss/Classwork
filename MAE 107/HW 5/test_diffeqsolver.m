function test_diffeqsolver

%
% Objective:
%   Solve a diff eq using a tridiagonal matrix
% input variables:
%   none
%
% output variables:
%   none
%
% functions called:
%   diffeqsolver - computes x for a tridiagonal augmented matrix
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
% Define n values
%
n = 2.^(2:5);
%
% find x and u values
%
[xn4, un4] = diffeqsolver(n(1),2,3);
[xn8, un8] = diffeqsolver(n(2),2,3);
[xn16, un16] = diffeqsolver(n(3),2,3);
[xn32, un32] = diffeqsolver(n(4),2,3);
%
% Plot xn and un for each iteration
%
% Create figure
%
figure(1); hold on;
%
% Plot n=4
%
plot(xn4,un4);
%
% Plot n=8
%
plot(xn8,un8);
%
% Plot n=16
%
plot(xn16,un16);
%
% Plot n=32
%
plot(xn32,un32);
%
% Lables
%
xlabel('x');
ylabel('u');
legend('n = 4', 'n = 8','n = 16','n = 32');
title('application of tridiagonal matrix to diff eqs');