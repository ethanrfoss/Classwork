%% Asignment 1.3
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

%for i=1:(301-150)
    Ui = squeeze(mean(u(151:end,:,:)));
    Vi = squeeze(mean(v(151:end,:,:)));
    
    
    pcolor(x,y,Ui);
    cu = colorbar;
    %cu.Limits = [min(min(min(u))) max(max(max(u)))];
    shading interp;
    rectangle('Position',[-0.5 -0.5 1 1],'Curvature',[1 1],'LineStyle','none','FaceColor',[1 1 1]);
    
    %[X,Y] = meshgrid(x(:,1)',y(1,:));
    %quiver(x,y,Ui,Vi);
    streamline(x',y',Ui',Vi',-5*ones(15,1),linspace(-5,5,15));

    title('Time Averaged u(t) with streamlines');
    xlabel('x');
    ylabel('y');

%end

