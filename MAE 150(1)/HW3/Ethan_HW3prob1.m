clear all; clc;
%% 1a
%Degrees converted to radians for simplicity
syms f1(x) f2(x) f3(x) f4(x) f5(x) f6(x) rho(theta);
f1(x) = 5*(1-cos(pi*x/(pi/3)));
f2(x) = 10;
f3(x) = 15*((x-(5*pi/9))/(7*pi/18)-sin(2*pi*(x-(5*pi/9))/(7*pi/18))/(2*pi))+10;
f4(x) = 25;
f5(x) = -25*(10*(x-(10*pi/9))^3/(5*pi/9)^3-15*(x-(10*pi/9))^4/(5*pi/9)^4+6*(x-(10*pi/9)).^5/(5*pi/9)^5)+25;
f6(x) = 0;
f(x) = piecewise(0<x<pi/3,f1(x),pi/3<x<5*pi/9,f2(x),5*pi/9<x<17*pi/18,f3(x),17*pi/18<x<10*pi/9,f4(x),10*pi/9<x<15*pi/9,f5(x),15*pi/9<x<2*pi,f6(x));
v(x) = diff(f)*(1000*pi/60);
a(x) = diff(v)*(1000*pi/60);
j(x) = diff(a)*(1000*pi/60);
rho(theta) = f(theta);
figure(1); hold on;
title("1a");
subplot(2,2,1);
fplot(f(x),[0,2*pi]);
xlabel("theta(radians)");
ylabel("displacement(mm)");
subplot(2,2,2);
fplot(v(x),[0,2*pi]);
xlabel("theta(radians)");
ylabel("velocity(mm/s)");
subplot(2,2,3);
fplot(a(x),[0,2*pi]);
xlabel("theta(radians)");
ylabel("acceleration(mm/s^2)");
subplot(2,2,4);
fplot(j(x),[0,2*pi]);
xlabel("theta(radians)");
ylabel("jerk(mm/s^3)");

%% 1b
figure(2); hold on;
pangle = atand((v/(1000*pi/60))/(35+f));
fplot(pangle,[0,2*pi]);
xlabel('theta(radians)');
ylabel('pressure angle(degrees)');
title('Pressure Angle Plot');
    
%% 1c
figure(3); polaraxes; hold on;
theta = 0:.01:2*pi;
r = f(theta);
polarplot(theta, r+30);%cam contour
polarplot(theta, r+35);%pitch curve circle
polarplot(theta, 35*ones(1,length(theta)));%prime circle
polarplot(theta, 30*ones(1,length(theta)));%base circle
legend('Cam Contour','Pitch Curve','Prime Circle','base circle');
title('1c');

%% 1d
figure(4); polaraxes;hold on;
for i=1:9
    subplot(3,3,i);
    polarplot(theta+i*pi/6, r+30);%cam contour
    hold on;
    polarplot(theta+i*pi/6, 35*ones(1,length(theta)));%prime circle
    title(sprintf('theta = %d',i*30));
end
