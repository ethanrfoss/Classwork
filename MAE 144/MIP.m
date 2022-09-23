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
h = figure(1);

t = 0:.01:10;
theta = step(1*minreal(D2*P*D1*G/(1+D1*G+D1*G*G2*D2*P))*tf([1 0],[1 0 1]),0:.01:10);
phi = step(1*T2*tf([1 0],[1 0 1]),0:.01:10);
u = phi*r;
clf;
hold on;
rec = rectangle('Position',[u(1)-r 0 2*r 2*r],'Curvature',[1 1]);
pen = plot([u(1) 2*l*sin(theta(1))+u(1)],[r r+2*l*cos(theta(1))]);
philine = plot([u(1) u(1)+r*cos(phi(1))],[r r-r*sin(phi(1))]);
%vert = plot([u(1) u(1)],[0 0.4080],'--');
leg = legend([pen philine],sprintf('theta = %.4f deg',theta(1)*180/pi),sprintf('phi = %.4f deg',phi(1)*180/pi));
hold off;
axis([-0.2040 0.2040 0 0.4080]);
axis equal;
tit = title(sprintf('Pendulum Animation t = 0'));
for i = 1:length(phi)
    set(rec,'Position',[u(i)-r 0 2*r 2*r]);
    set(pen,'XData',[u(i) 2*l*sin(theta(i))+u(i)],'YData',[r r+2*l*cos(theta(i))]);
    set(philine,'XData',[u(i) u(i)+r*cos(phi(i))],'YData',[r r-r*sin(phi(i))]);
    set(leg,'String', { sprintf('theta = %.4f deg',theta(i)*180/pi)   sprintf('phi = %.4f deg',phi(i)*180/pi)});
    set(tit,'String',sprintf('Pendulum Animation t = %.2f',t(i)));
    %set(vert,'XData', [u(i) u(i)], 'YData', [0 0.4080]);
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1 
          imwrite(imind,cm,'PendulumAnimation.gif','gif','DelayTime',0.01, 'Loopcount',inf); 
      else 
          imwrite(imind,cm,'PendulumAnimation.gif','gif','DelayTime',0.01,'WriteMode','append'); 
      end
end