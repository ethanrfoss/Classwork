clear all; clc;
%Degrees converted to radians for simplicity
%% 2a
syms f1(x) f2(x) f3(x) f4(x);
f1(x) = 22*(x/(5*pi/9)-sin(2*pi*x/(5*pi/9))/(2*pi));
f2(x) = 22;
f3(x) = 22-22*(10*(x-(11*pi/12))^3/(23*pi/36)^3 - 15*(x-(11*pi/12))^4/(23*pi/36)^4 + 6*(x-(11*pi/12))^5/(23*pi/36)^5);
f4(x) = 0;
f(x) = piecewise(0<=x<5*pi/9,f1(x),5*pi/9<=x<11*pi/12,f2(x),11*pi/12<=x<14*pi/9,f3(x),14*pi/9<=x<2*pi,f4(x));
v(x) = diff(f)*(50*pi);
a(x) = diff(v)*(50*pi);
theta = 0:pi/180:2*pi;
dis = f(theta);
vel = v(theta);
acc = a(theta);

figure(1); hold on;
title("2a");
subplot(1,3,1);
plot(theta,dis);
xlabel("theta(degrees)");
ylabel("displacement(mm)");
subplot(1,3,2);
plot(theta,vel);
xlabel("theta(degrees)");
ylabel("velocity(mm/s)");
subplot(1,3,3);
plot(theta,acc);
xlabel("theta(degrees)");
ylabel("acceleration(mm/s^2)");

%% 2b
figure(2); hold on;
pangle = atand((v/(50*pi))/(68+f));
fplot(pangle,[0,2*pi]);

%% 2c

figure(3); hold on;
[xcam,ycam] = pol2cart(theta+pi/2,(dis+60));
plot(xcam,ycam); %cam contour
[xpitch,ypitch] = pol2cart(theta+pi/2,(dis+68));
plot(xpitch,ypitch); %pitch curve
[xbase,ybase] = pol2cart(theta, 60*ones(1,length(theta)));
plot(xbase,ybase); %base circle
[xprime,yprime] = pol2cart(theta, 68*ones(1,length(theta)));
plot(xprime,yprime); %prime circle
legend('Cam contour','Pitch Curve','Base Circle','Prime Circle');
axis([-100 100 -100 100]);

%% 2d
maxvel = max(vel);

%% 2e
maxacc = max(acc);