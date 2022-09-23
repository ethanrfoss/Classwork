clear all; close all; clc;

%% 4.1
xvals = (0:.1:30)/100;
deltax = @(x) 5*sqrt(1.5116*10^-5*x/.5);
figure(1);
plot(xvals,deltax(xvals));
title('Boundary Layer Thickness vs. Length')
xlabel('Distnace(m)');
ylabel('Boundary Layer Thickness(m)');
