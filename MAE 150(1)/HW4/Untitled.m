t2 = 0:.1:2*pi;
t4 = zeros(1,length(t2));
error = .1;
f1 = @(theta2,theta4) 466+540*cos(theta4)-216*cos(theta2)-180*cos(theta4)*cos(theta2)-150*sin(theta4)-180*sin(theta4)*sin(theta2)+60*sin(theta2);
df1 = @(theta2,theta4) -540*sin(theta4)+180*sin(theta4)*cos(theta2)-150*cos(theta4)-180*cos(theta4)*sin(theta2);
for i = 1:length(t2)
    t4temp = zeros(1,100);
    t4temp(1) = pi*6/5;
    j = 1;
    h = 1;
    while (abs(h)>=error)
        h = f1(t2(i),t4temp(j))/df1(t2(i),t4temp(j));
        t4temp(j+1) = t4temp(j)-h;
    end
    t4(i) = t4temp(j);
end
figure(1);
plot(t2,t4);
xlabel('Crank angle');
ylabel('Angle of Rocker');
title('2a - Input Crank Angle vs Rocker Angle');