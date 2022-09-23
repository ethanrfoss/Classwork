%% Asignment 2.3
clear all; close all; clc;

%% Open .mat file
load('cylinder_Re100.mat')
% dt = .0050
% u = 301x200x100
% v = 301x200x100
% x = 200x100
% y = 200x100

%% Figure
dx = .1005;
dy = .1010;
figure(1); hold on;

    
for i=1:301
    Up=squeeze(u(i,:,:));
    Vp=squeeze(v(i,:,:));

    
    w = ddx_central(Vp,dx)-ddy_central(Up,dy);
    pcolor(x,y,w);
    cu = colorbar;
    %cu.Limits = [min(min(min(u))) max(max(max(u)))];
    shading interp;
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);

    title('w');
    xlabel('x');
    ylabel('y');
    drawnow;
end
