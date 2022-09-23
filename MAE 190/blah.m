%% Asignment 1.5
clear all; close all; clc;

%% Open .mat file
load('cylinder_Re100.mat')
% dt = .0050
% u = 301x200x100
% v = 301x200x100
% x = 200x100
% y = 200x100

%% Figure

figure(1); hold on;

Ubar = squeeze(mean(u(151:301,:,:)));
Vbar = squeeze(mean(v(151:301,:,:)));
    
for i=151:301
    Usquared(i-150,:,:) = (squeeze(u(i,:,:))-Ubar).^2;
    Vsquared(i-150,:,:) = (squeeze(v(i,:,:))-Vbar).^2;
end
    
    k = .5*(squeeze(mean(Usquared))+squeeze(mean(Vsquared)));
    pcolor(x,y,k);
    cu = colorbar;
    %cu.Limits = [min(min(min(u))) max(max(max(u)))];
    shading interp;
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);

    title('k');
    xlabel('x');
    ylabel('y');
    drawnow;

