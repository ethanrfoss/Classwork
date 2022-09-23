%% Assignment 4_3
clear all; close all; clc;
%% Open .mat file
load('MAE190_4_3','u','v','xx','yy');
% u 200x100
% v 200x100
% xx 200x100
% yy 200x100

%% Constants
dx = .0101; dy = .0101;
d = .0005;
%% Calculate

Un = u;
Unbar = u;
Vn = v;
Vnbar = v;
n = 1;
t = 0;
figure(1);
subplot(2,1,1);
pu = pcolor(xx,yy,Un); shading interp; colorbar;
subplot(2,1,2);
pv = pcolor(xx,yy,Vn); shading interp; colorbar;
while t<=2
    
    dt = 1/(max(max(Un))/dx+max(max(Vn))/dy);
    t = t + dt;
    
    Unbar = Un-Un.*dt.*ddx_fwd_periodic(Un,dx)-Vn.*dt.*ddy_fwd(Un,dy) + d*dt*(d2dx2periodic(Un,dx)+d2dy2(Un,dy));
    Unbar(:,1) = 0; Unbar(:,end) = 0;
    Un = .5*((Un+Unbar)-Un.*dt.*ddx_bwd_periodic(Unbar,dx)-Vn.*dt.*ddy_bwd(Unbar,dy) + d*dt*(d2dx2periodic(Un,dx)+d2dy2(Un,dy)));
    
    Vnbar = Vn-Un.*dt.*ddx_fwd_periodic(Vn,dx)-Vn.*dt.*ddy_fwd(Vn,dy) + d*dt*(d2dx2periodic(Vn,dx)+d2dy2(Vn,dy));
    Vnbar(:,1) = 0; Vnbar(:,end) = 0;
    Vn = .5*((Vn+Vnbar)-Un.*dt.*ddx_bwd_periodic(Vnbar,dx)-Vn.*dt.*ddy_bwd(Vnbar,dy) + d*dt*(d2dx2periodic(Vn,dx)+d2dy2(Vn,dy)));
    
    n = n+1;
    
    subplot(2,1,1);
    set(pu,'CData',Un); shading interp; colorbar;
    drawnow;
    title(sprintf('2D Baterman-Burgers Flow in u direction(t = %d)',t));
    
    subplot(2,1,2);
    set(pv,'CData',Vn); shading interp; colorbar;
    drawnow;
    title(sprintf('2D Baterman-Burgers Flow in v direction(t = %d)',t));
end
