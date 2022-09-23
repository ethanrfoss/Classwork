
%% Water
rho = 996;
v = .87*10^-6;
k = .611;
cp = 4178;
D = .01;

V = .01:.01:100;
hc = zeros(1, length(V));
Re = zeros(1, length(V));
for i = 1:length(V)
    Re(i) = D*V(i)/v;
    if Re(i)<10^4
        hc(i) = 3.66*k/D;
    else
        hc(i) = .023*V(i)^.8*k^.6*(rho*cp)^.4/(D^.2*v^.4);
    end
end

figure(1); hold on;
plot(V,hc);
%plot(V,Re);
title('Heat Transfer Coefficient of Water over a Range of Fluid Velocities');
xlabel('V(m/s)');
ylabel('hc(W/m^2K)');


%% Air
rho = 1.177;
v = 15.7*10^-6;
k = .0267;
cp = 1005;
D = .01;

V = .01:.01:100;
hc = zeros(1, length(V));
Re = zeros(1, length(V));
for i = 1:length(V)
    Re(i) = D*V(i)/v;
    if Re(i)<10^4
        hc(i) = 3.66*k/D;
    else
        hc(i) = .023*V(i)^.8*k^.6*(rho*cp)^.4/(D^.2*v^.4);
    end
end

figure(2); hold on;
plot(V,hc);
%plot(V,Re);
title('Heat Transfer Coefficient of Air over a Range of Fluid Velocities');
xlabel('V(m/s)');
ylabel('hc(W/m^2K)');

%% 2.2

x = 0:.01:1;
hc = zeros(1,length(x));
for i =1:length(x)
    if x(i)<.71422
        hc(i) = 1.07*(20/x(i))^.25;
    else
        hc(i) = 1.3*(20)^.333333333;
    end
end

figure(3); hold on;
plot(hc,x);
%plot(V,Re);
title('Heat Transfer Coefficient of Air over a 1m Vertical Wall');
xlabel('hc(W/m^2K)');
ylabel('x(m)');
 