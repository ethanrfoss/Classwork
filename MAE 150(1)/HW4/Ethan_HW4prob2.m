clear all; clc;

%% 2a
t2 = 0:.05:2*pi;
t4 = zeros(1,length(t2));
error = .05;
f1 = @(theta2,theta4) 466+540*cos(theta4)-216*cos(theta2)-180*cos(theta4)*cos(theta2)-150*sin(theta4)-180*sin(theta4)*sin(theta2)+60*sin(theta2);
df1 = @(theta2,theta4) -540*sin(theta4)+180*sin(theta4)*cos(theta2)-150*cos(theta4)-180*cos(theta4)*sin(theta2);
for i = 1:length(t2)
    t4temp = zeros(1,100);
    t4temp(1) = pi*3/5;
    j = 1;
    h = 1;
    while (abs(h)>=error)
        h = f1(t2(i),t4temp(j))/df1(t2(i),t4temp(j));
        t4temp(j+1) = t4temp(j)-h;
        j = j+1;
    end
    t4(i) = t4temp(j);
end
figure(1);
plot(t2*180/pi,t4*180/pi);
xlabel('Crank angle(degrees)');
ylabel('Angle of Rocker(degrees)');
title('2a - Input Crank Angle vs Rocker Angle');

%% 2b
t3 = acos((18+15*cos(t4)-6*cos(t2))/12);
xp = 6*cos(t2)-5*sin(t3);
yp = 6*sin(t2)+5*cos(t3);
figure(2); hold on;
plot(xp,yp);
xlabel('x');
ylabel('y');
title('2b - Path of point P');

%% 2c
xq = 18+22*cos(t4);
yq = -5+22*sin(t4);
d = sqrt((xq-xp).^2+(yq-yp).^2);
figure(3); hold on;
title('2c - Crank angle vs. Distance between P and Q');
xlabel('Crank angle(degrees)');
ylabel('Distance between P and Q(cm)');
plot(t2*180/pi,d);

%% 2d
w2 = 100*2*pi/60;
w4 = (6*w2*sin(t2)./tan(t3)-6*w2*cos(t2))./(15*sin(t4)./tan(t3)-15*cos(t4));
figure(4); hold on;
title('2d - Crank angle vs. Angular velocity of long link');
xlabel('Crank angle(degrees)');
ylabel('Angular velocity of long link(rad/s)');
plot(t2*180/pi,w4);

%% 2e

w3 = (15*w4.*sin(t4)-6*w2*sin(t2))./(12*sin(t3));
vp = sqrt((-6*w2*sin(t2)-5*w3.*cos(t3)).^2+(6*w2*cos(t2)-5*w3.*sin(t3)).^2);
figure(5); hold on;
title('2e - Crank angle vs. Translational velocity of point P');
xlabel('Crank angle(degrees)');
ylabel('Translational Velocity of P(cm/s)');
plot(t2*180/pi,vp);

%% Moving CAM plot

figure(6); hold on;
axis([-20 20 -20 20]);
plot(xp,yp,'--g');
plot(xq,yq,'--g');
plot(0,0,'or');
xl23 = 6*cos(t2);
yl23 = 6*sin(t2);
xl34 = 18+15*cos(t4);
yl34 = -5 + 15*sin(t4);
p1 = plot([0 xl23(1)],[0 yl23(1)],'r'); %l2
p2 =plot([xl23(1) xl34(i)],[yl23(1) yl34(1)],'r'); %l3
p3 = plot([18 xq(1)],[-5 yq(1)],'r');%l4
p4 = plot([xl34(1) xp(1)],[yl34(1) yp(1)],'r'); %p
for i = 1:length(t2)
    delete(p1);
    delete(p2);
    delete(p3);
    delete(p4);
    p1 = plot([0 xl23(i)],[0 yl23(i)],'-r'); %l2
    p2 =plot([xl23(i) xl34(i)],[yl23(i) yl34(i)],'-r'); %l3
    p3 = plot([18 xq(i)],[-5 yq(i)],'-r');%l4
    p4 = plot([xl34(i) xp(i)],[yl34(i) yp(i)],'-r'); %p
    plot(xl23(i),yl23(i),'-b');
    pause(.1);
end