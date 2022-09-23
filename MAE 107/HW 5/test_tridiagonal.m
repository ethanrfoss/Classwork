function test_tridiagonal

%
% Objective:
%   Solve a tridiagonal matrix
% input variables:
%   none
%
% output variables:
%   none
%
% functions called:
%   tridiagonal - computes x for a tridiagonal augmented matrix
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
% Define matrix and determine x
%
u = [1 3 3 -1 -3];
d = [2 6 6 8 9 8];
l = [2 1 2 -2 1];
b = [1 2 3 4 5 6];
x = tridiagonal(u,d,l,b)
A = [ 2 1 0 0 0 0; 2 6 3 0 0 0; 0 1 6 3 0 0; 0 0 2 8 -1 0; 0 0 0 -2 9 -3; 0 0 0 0 1 8];
%
% Check accuracy by computing euclidian norm
%
A*x-b'