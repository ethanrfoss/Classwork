%% Asignment 3.2 part 2
clear all; close all; clc;
%% Open .mat file
load('mae190','T','xx','yy');
% T = 400x200
% xx = 400x200
% yy = 400x200

%% Figure
figure(1); hold on;
c = [1 1];
dtmax = .0025;
dt = .001; dx = .005; dy = .005;
t = 0:dt:2;

figure(1);
Tn = zeros(400,200,length(t));
Tn(:,:,1) = T;
p = pcolor(xx,yy,Tn(:,:,1)); shading interp; colorbar;
for n = 1:length(t)-1
    
    Tn(:,:,n+1) = Tn(:,:,n) -c(1)*dt*ddx_central(Tn(:,:,n),dx)-c(2)*dt*ddy_central(Tn(:,:,n),dy);
    
    set(p,'CData',Tn(:,:,n+1)); shading interp; colorbar;
    drawnow;
    title(sprintf('2D Heat Equation(t = %d)',t(n+1)));
end

