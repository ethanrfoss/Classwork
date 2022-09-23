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
dx = .0042; dy = .0034;
%% Calculate

mu = sutherland(T);

txx1 = 2*mu.*(ddx_fwd(u,dx)-1/3*(ddx_fwd(u,dx)+ddy_central(u,dy)));
tyy1 = 2*mu.*(ddy_central(v,dy)-1/3*(ddx_fwd(v,dx)+ddy_central(v,dy)));
txy1 = mu.*(ddy_central(u,dy)+ddx_fwd(v,dx));

txx2 = 2*mu.*(ddx_central(u,dx)-1/3*(ddx_central(u,dx)+ddy_bwd(u,dy)));
tyy2 = 2*mu.*(ddy_bwd(v,dy)-1/3*(ddx_central(v,dx)+ddy_bwd(v,dy)));
txy2 = mu.*(ddy_bwd(u,dy)+ddx_central(v,dx));

%% Plot

figure(1); hold on;

subplot(2,3,1);
pcolor(xx,yy,txx1); shading interp;
cb = colorbar; ylabel(cb,'\tauxx [N/m^2]'); caxis([-.5 .5]);
title('\tauxx(fwd x, cen y)'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(2,3,2);
pcolor(xx,yy,tyy1); shading interp;
cb = colorbar; ylabel(cb,'\tauyy [N/m^2]'); caxis([-.5 .5]);
title('\tauyy(fwd x, cen y)'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(2,3,3);
pcolor(xx,yy,txy1); shading interp;
cb = colorbar; ylabel(cb,'\tauxy [N/m^2]'); caxis([-.5 .5]);
title('\tauxy(fwd x, cen y)'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(2,3,4);
pcolor(xx,yy,txx2); shading interp;
cb = colorbar; ylabel(cb,'\tauxx [N/m^2]'); caxis([-.5 .5]);
title('\tauxx(cen x, bwd y)'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(2,3,5);
pcolor(xx,yy,tyy2); shading interp;
cb = colorbar; ylabel(cb,'\tauyy [N/m^2]'); caxis([-.5 .5]);
title('\tauyy(cen x, bwd y)'); xlabel('x'); ylabel('y');
axis equal tight;

subplot(2,3,6);
pcolor(xx,yy,txy2); shading interp;
cb = colorbar; ylabel(cb,'\tauxy [N/m^2]'); caxis([-.5 .5]);
title('\tauxy(cen x, bwd y)'); xlabel('x'); ylabel('y');
axis equal tight;

