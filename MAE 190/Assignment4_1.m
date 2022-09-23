%% Assignment 4_1
clear all; close all; clc;
%% Open .mat file
load('supersonicJetLES_xyPlane','rho','u','v','T','xx','yy');
% rho 300x75
% T 300x75
% u 300x75
% v 300x75
% xx 300x75
% yy 300x75
%% Constants
cp = 1005; %J/kgk
cv = 718; %J/kgk
R = cp-cv; %J/kgk

%% Calculate

U = prim2cons(rho,u,v,T,cv);
[rho,u,v,T,p,e,Et] = cons2prim(U,R,cv);
mu = sutherland(T);
%% Plot

figure(1); hold on;

subplot(4,2,1);
pcolor(xx,yy,rho); shading interp;
cb = colorbar; ylabel(cb,'\rho [kg/m^3]')
title('\rho'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(4,2,2);
pcolor(xx,yy,squeeze(U(2,:,:))); shading interp;
cb = colorbar; ylabel(cb,'\rhou [kg/s/m^2]')
title('\rhou'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(4,2,3);
pcolor(xx,yy,squeeze(U(3,:,:))); shading interp;
cb = colorbar; ylabel(cb,'\rhov [kg/s/m^2]')
title('\rhov'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(4,2,4);
pcolor(xx,yy,T); shading interp;
cb = colorbar; ylabel(cb,'T [K]')
title('T'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(4,2,5);
pcolor(xx,yy,p); shading interp;
cb = colorbar; ylabel(cb,'p [N/m^2]')
title('p'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(4,2,6);
pcolor(xx,yy,e); shading interp;
cb = colorbar; ylabel(cb,'e [Nm/kg]')
title('e'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(4,2,7);
pcolor(xx,yy,Et); shading interp;
cb = colorbar; ylabel(cb,'Et [Nm]')
title('Et'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(4,2,8);
pcolor(xx,yy,mu); shading interp;
cb = colorbar; ylabel(cb,'\mu [Ns/m^2]')
title('\mu'); xlabel('x'); ylabel('y');
axis equal tight;
