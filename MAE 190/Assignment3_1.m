%% Asignment 3.1
clear all; close all; clc;
%% Open .mat file
load('mae190','T','xx','yy');
% T = 400x200
% xx = 400x200
% yy = 400x200

%% Figure
figure(1); hold on;
dtmax = 3.125*10^-6;
dt = 1*10^-6; dx = .005; dy = .005;
a = 2;
t = 0:dt:.001;

Tn = zeros(400,200,length(t));
Tn(:,:,1) = T;
for n = 1:length(t)-1
    pcolor(xx,yy,Tn(:,:,n)); shading interp; colorbar;
    
    Tn(:,:,n+1) = Tn(:,:,n) + a*dt*d2dx2periodic(Tn(:,:,n),dx) + a*dt*d2dy2periodic(Tn(:,:,n),dy);
    
    drawnow;
    title(sprintf('2D Heat Equation(t = %d)',t(n)));
end

