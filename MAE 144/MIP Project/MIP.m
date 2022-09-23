%% MIP Project

%% Paramters of MIP

st = .003; %[Nm]
G = 35.555555555555555; 
Im = 3.6*10^-8; %[Kgm^2]
r = 34/1000; %[mm]
mw = 27/1000; %[kg]
mb = 180/1000; %[kg]
l = 47.7/1000; %[m]
Ib = 2.63*10^-4; %[kgm^2]
Iw = 2*(mw*r^2/2+G^2*Im); %[kgm^2]
g = 9.81; %[m/s^2]
k = st/1760; %[Nms]

%% Plant Constants

alpha1 = (Iw + (mw+mb)*r^2)*(2*G*st)+(2*G*st)*(mb*r*l);
beta1 = -(Iw+(mw+mb)*r^2)*(Ib+mb*l^2)+(mb*r*l)^2;
beta2 = (2*k*G^2)*(-(Iw+(mw+mb)*r^2)-(Ib+mb*l^2)-2*(mb*r*l));
beta3 = (Iw +(mw+mb)*r^2)*(mb*g*l);
beta4 = (2*k*G^2)*(mb*g*l);

%% Plant

G = tf([alpha1 0],[beta1 beta2 beta3 beta4]);
minreal(G);
[b a] = tfdata(G, 'v');
b = tf(b,1);
a = tf(a,1);


%% D1 (found using sisotool)

D1 = -100*tf([.1806 1],[6.7 1]);
y = -100*tf([.18 1],1);
x = tf([6.7 1],1);
wc = 38.3; %[rad/s]
h = 1/100; %[s]
f = 2*(1-cos(wc*h))/(wc*h*sin(wc*h));
syms z;
D1z = -100*((2*.1806+f*h)*z+(f*h-2*.1806))/((2*6.7+f*h)*z+(f*h-2*6.7));

%% P (from step resonse)

P = 1/1.21;
%% T1

T1 = b*y/(a*x+b*y);
%% Impulse Response

impulse(T)

%% G2

alpha12 = -(Ib+mb*l^2)-(mb*r*l);
alpha22 = 0;
alpha32 = (mb*g*l);
beta12 = (Iw + (mw+mb)*r^2)+(mb*r*l);

G2 = tf([alpha12 0 alpha32],[beta12 0 0]);

%% D2 (found using sisotool)

D2 = .0071465*tf([1 1],[.1 1]);
wc = 1.21; %[rad/s]
h = 1/20; %[s]
f = 2*(1-cos(wc*h))/(wc*h*sin(wc*h));
syms z;
D2z = .0071465*((2+f*h)*z+(f*h-2))/((.2+f*h)*z+(f*h-.2));

%% T2

GD2 = D2*P*T1*G2;

[num dem] = tfdata(GD2, 'v');
num = tf(num,1);
dem = tf(dem,1);

T2 = num/(dem+num);

%% Simulation

figure(1);

theta = step(*D2*P/(1+D1*D1+D2*P),0:.01:10);
phi = step(*T2,0:.01:10);
u = phi*r;
clf;
hold on;
rec = rectangle('Position',[u(1)-r 0 2*r 2*r],'Curvature',[1 1]);
pen = plot([u(1) 2*l*sin(theta(1))+u(1)],[r r+2*l*cos(theta(1))]);
hold off;
axis([-0.2040 0.2040 0 0.4080]);
axis equal;
for i = 1:length(phi)
    set(rec,'Position',[u(i)-r 0 2*r 2*r]);
    set(pen,'XData',[u(i) 2*l*sin(theta(i))+u(i)],'YData',[r r+2*l*cos(theta(i))]);
    pause(.005)
end