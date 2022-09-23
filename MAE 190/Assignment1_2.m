%% Asignment 1.2
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



up = subplot(2,1,1);
vp = subplot(2,1,2);
%for i=1:(301-150)
    up = subplot(2,1,1);
    
    Ui = squeeze(mean(u(151:301,:,:)));
    pcolor(x,y,Ui);
    cu = colorbar;
    %cu.Limits = [min(min(min(u))) max(max(max(u)))];
    shading interp;
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);

    title('ubar(transient discarded)');
    xlabel('x');
    ylabel('y');

    vp = subplot(2,1,2);

    Vi = squeeze(mean(v(151:301,:,:)));
    pcolor(x,y,Vi);
    cv = colorbar;
    %cv.Limits = [min(min(min(v))) max(max(max(v)))];
    shading interp;
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);

    title('vbar(transient discarded)');
    xlabel('x');
    ylabel('y');
    drawnow;
%end

