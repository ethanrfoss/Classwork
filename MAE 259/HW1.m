
%% Homework 1 MAE 259
% Ethan Foss

%% Problem 2

syms s L; % L = h*lambda

P = [1, -(1+23/12*L),  16/12*L, -5/12*L];

r = roots(P);
r1 = r(1);

sympref('PolynomialDisplayStyle','ascend');
disp("Polynomial of First Root: ")
series(r1,L,'Order',6)

r2 = r(2);
disp("First Spurious Root: ")
series(r2,L,'Order',6)

r3 = r(3);
disp("Second Spurious Root: ")
series(r3,L,'Order',6)