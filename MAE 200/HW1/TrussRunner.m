StudentIDnumber=16129635; configuration=mod(StudentIDnumber,4)+1

%% 4.1

g = [1 0 1 1 2 0 2 1];
xa = 0; ya = 0;
xb = 0; yb = 1;
xg = 3; yg = 1;

f = TrussForces(g)

%% 4.2

l = [sqrt((g(2)-ya)^2+(g(1)-xa)^2);sqrt((g(4)-ya)^2+(g(3)-xa)^2);sqrt((g(4)-yb)^2+(g(3)-xb)^2);sqrt((g(4)-g(2))^2+(g(3)-g(1))^2);...
    sqrt((g(6)-g(2))^2+(g(5)-g(1))^2); sqrt((g(8)-g(2))^2+(g(7)-g(1))^2); sqrt((g(8)-g(4))^2+(g(7)-g(3))^2); sqrt((g(8)-g(6))^2+(g(7)-g(5))^2);...
    sqrt((yg-g(6))^2+(xg-g(5))^2); sqrt((yg-g(8))^2+(xg-g(7))^2)];

rho = 5;
a = -9.81;
W = -1000;
z = [0; .5*(rho*l(1)+rho*l(5)+rho*l(6))*a; 0; .5*(rho*l(2))*a; 0; .5*(rho*l(5)+rho*l(9))*a; 0; .5*(rho*l(6))*a; 0; .5*(rho*l(9))*a+W];

m = [.5*(rho*l(1)+rho*l(5)+rho*l(6)); .5*(rho*l(1)+rho*l(5)+rho*l(6)); .5*(rho*l(2)); .5*(rho*l(2)); .5*(rho*l(5)+rho*l(9)); .5*(rho*l(5)+rho*l(9)); .5*(rho*l(6)); .5*(rho*l(6)); .5*(rho*l(9)); .5*(rho*l(9))];


[xe,fe] = TrussLoading(g,z);

Scale = 100;
figure(1); hold on;
plot([xa g(1)],[ya g(2)],'g');
plot([xa g(3)],[ya g(4)],'g');
plot([xb g(3)],[yb g(4)],'r');
plot([g(1) g(3)],[g(2) g(4)],'r');
plot([g(1) g(5)],[g(2) g(6)],'g');
plot([g(1) g(7)],[g(2) g(8)],'g');
plot([g(3) g(7)],[g(4) g(8)],'r');
plot([g(5) g(7)],[g(6) g(8)],'r');
plot([g(5) xg],[g(6) yg],'g');
plot([g(7) xg],[g(8) yg],'r');
axis equal; axis([-.5 3.5 -.5 1.5]);

plot([xa g(1)+Scale*xe(1)],[ya g(2)+Scale*xe(2)],'g');
plot([xa g(3)+Scale*xe(3)],[ya g(4)+Scale*xe(4)],'g');
plot([xb g(3)+Scale*xe(3)],[yb g(4)+Scale*xe(4)],'r');
plot([g(1)+Scale*xe(1) g(3)+Scale*xe(3)],[g(2)+Scale*xe(2) g(4)+Scale*xe(4)],'r');
plot([g(1)+Scale*xe(1) g(5)+Scale*xe(5)],[g(2)+Scale*xe(2) g(6)+Scale*xe(6)],'g');
plot([g(1)+Scale*xe(1) g(7)+Scale*xe(7)],[g(2)+Scale*xe(2) g(8)+Scale*xe(8)],'g');
plot([g(3)+Scale*xe(3) g(7)+Scale*xe(7)],[g(4)+Scale*xe(4) g(8)+Scale*xe(8)],'r');
plot([g(5)+Scale*xe(5) g(7)+Scale*xe(7)],[g(6)+Scale*xe(6) g(8)+Scale*xe(8)],'r');
plot([g(5)+Scale*xe(5) xg+Scale*xe(9)],[g(6)+Scale*xe(6) yg+Scale*xe(10)],'g');
plot([g(7)+Scale*xe(7) xg+Scale*xe(9)],[g(8)+Scale*xe(8) yg+Scale*xe(10)],'r');

title('Original and Deformed Truss(100 Magnification');

%% 4.3

[K C D] = TrussK(g);

M = diag(m);
L = chol(M);

A = L^-1*K*(L')^-1;

[S,E] = eig(A);

Y = (L')^-1*S;
[YM,EM] = eig(M^-1*K);



figure(2); hold on;
Scale = 1;
plot([xa g(1)+Scale*Y(1,6)],[ya g(2)+Scale*Y(2,6)],'r');
plot([xa g(3)+Scale*Y(3,6)],[ya g(4)+Scale*Y(4,6)],'r');
plot([xb g(3)+Scale*Y(3,6)],[yb g(4)+Scale*Y(4,6)],'r');
plot([g(1)+Scale*Y(1,6) g(3)+Scale*Y(3,6)],[g(2)+Scale*Y(2,6) g(4)+Scale*Y(4,6)],'r');
plot([g(1)+Scale*Y(1,6) g(5)+Scale*Y(5,6)],[g(2)+Scale*Y(2,6) g(6)+Scale*Y(6,6)],'r');
plot([g(1)+Scale*Y(1,6) g(7)+Scale*Y(7,6)],[g(2)+Scale*Y(2,6) g(8)+Scale*Y(8,6)],'r');
plot([g(3)+Scale*Y(3,6) g(7)+Scale*Y(7,6)],[g(4)+Scale*Y(4,6) g(8)+Scale*Y(8,6)],'r');
plot([g(5)+Scale*Y(5,6) g(7)+Scale*Y(7,6)],[g(6)+Scale*Y(6,6) g(8)+Scale*Y(8,6)],'r');
plot([g(5)+Scale*Y(5,6) xg+Scale*Y(9,6)],[g(6)+Scale*Y(6,6) yg+Scale*Y(10,6)],'r');
plot([g(7)+Scale*Y(7,6) xg+Scale*Y(9,6)],[g(8)+Scale*Y(8,6) yg+Scale*Y(10,6)],'r');
axis equal; axis([-.5 3.5 -.5 1.5]);
title(sprintf('w = %f rad/s',sqrt(E(6,6))));


figure(3); hold on;
Scale = 1;
plot([xa g(1)+Scale*Y(1,7)],[ya g(2)+Scale*Y(2,7)],'r');
plot([xa g(3)+Scale*Y(3,7)],[ya g(4)+Scale*Y(4,7)],'r');
plot([xb g(3)+Scale*Y(3,7)],[yb g(4)+Scale*Y(4,7)],'r');
plot([g(1)+Scale*Y(1,7) g(3)+Scale*Y(3,7)],[g(2)+Scale*Y(2,7) g(4)+Scale*Y(4,7)],'r');
plot([g(1)+Scale*Y(1,7) g(5)+Scale*Y(5,7)],[g(2)+Scale*Y(2,7) g(6)+Scale*Y(6,7)],'r');
plot([g(1)+Scale*Y(1,7) g(7)+Scale*Y(7,7)],[g(2)+Scale*Y(2,7) g(8)+Scale*Y(8,7)],'r');
plot([g(3)+Scale*Y(3,7) g(7)+Scale*Y(7,7)],[g(4)+Scale*Y(4,7) g(8)+Scale*Y(8,7)],'r');
plot([g(5)+Scale*Y(5,7) g(7)+Scale*Y(7,7)],[g(6)+Scale*Y(6,7) g(8)+Scale*Y(8,7)],'r');
plot([g(5)+Scale*Y(5,7) xg+Scale*Y(9,7)],[g(6)+Scale*Y(6,7) yg+Scale*Y(10,7)],'r');
plot([g(7)+Scale*Y(7,7) xg+Scale*Y(9,7)],[g(8)+Scale*Y(8,7) yg+Scale*Y(10,7)],'r');
axis equal; axis([-.5 3.5 -.5 1.5]);
title(sprintf('w = %f rad/s',sqrt(E(7,7))));

figure(4); hold on;
Scale = 1;
plot([xa g(1)+Scale*Y(1,8)],[ya g(2)+Scale*Y(2,8)],'r');
plot([xa g(3)+Scale*Y(3,8)],[ya g(4)+Scale*Y(4,8)],'r');
plot([xb g(3)+Scale*Y(3,8)],[yb g(4)+Scale*Y(4,8)],'r');
plot([g(1)+Scale*Y(1,8) g(3)+Scale*Y(3,8)],[g(2)+Scale*Y(2,8) g(4)+Scale*Y(4,8)],'r');
plot([g(1)+Scale*Y(1,8) g(5)+Scale*Y(5,8)],[g(2)+Scale*Y(2,8) g(6)+Scale*Y(6,8)],'r');
plot([g(1)+Scale*Y(1,8) g(7)+Scale*Y(7,8)],[g(2)+Scale*Y(2,8) g(8)+Scale*Y(8,8)],'r');
plot([g(3)+Scale*Y(3,8) g(7)+Scale*Y(7,8)],[g(4)+Scale*Y(4,8) g(8)+Scale*Y(8,8)],'r');
plot([g(5)+Scale*Y(5,8) g(7)+Scale*Y(7,8)],[g(6)+Scale*Y(6,8) g(8)+Scale*Y(8,8)],'r');
plot([g(5)+Scale*Y(5,8) xg+Scale*Y(9,8)],[g(6)+Scale*Y(6,8) yg+Scale*Y(10,8)],'r');
plot([g(7)+Scale*Y(7,8) xg+Scale*Y(9,8)],[g(8)+Scale*Y(8,8) yg+Scale*Y(10,8)],'r');
axis equal; axis([-.5 3.5 -.5 1.5]);
title(sprintf('w = %f rad/s',sqrt(E(8,8))));

figure(5); hold on;
Scale = 1;
plot([xa g(1)+Scale*Y(1,9)],[ya g(2)+Scale*Y(2,9)],'r');
plot([xa g(3)+Scale*Y(3,9)],[ya g(4)+Scale*Y(4,9)],'r');
plot([xb g(3)+Scale*Y(3,9)],[yb g(4)+Scale*Y(4,9)],'r');
plot([g(1)+Scale*Y(1,9) g(3)+Scale*Y(3,9)],[g(2)+Scale*Y(2,9) g(4)+Scale*Y(4,9)],'r');
plot([g(1)+Scale*Y(1,9) g(5)+Scale*Y(5,9)],[g(2)+Scale*Y(2,9) g(6)+Scale*Y(6,9)],'r');
plot([g(1)+Scale*Y(1,9) g(7)+Scale*Y(7,9)],[g(2)+Scale*Y(2,9) g(8)+Scale*Y(8,9)],'r');
plot([g(3)+Scale*Y(3,9) g(7)+Scale*Y(7,9)],[g(4)+Scale*Y(4,9) g(8)+Scale*Y(8,9)],'r');
plot([g(5)+Scale*Y(5,9) g(7)+Scale*Y(7,9)],[g(6)+Scale*Y(6,9) g(8)+Scale*Y(8,9)],'r');
plot([g(5)+Scale*Y(5,9) xg+Scale*Y(9,9)],[g(6)+Scale*Y(6,9) yg+Scale*Y(10,9)],'r');
plot([g(7)+Scale*Y(7,9) xg+Scale*Y(9,9)],[g(8)+Scale*Y(8,9) yg+Scale*Y(10,9)],'r');
axis equal; axis([-.5 3.5 -.5 1.5]);
title(sprintf('w = %f rad/s',sqrt(E(9,9))));

figure(6); hold on;
Scale = 1;
plot([xa g(1)+Scale*Y(1,10)],[ya g(2)+Scale*Y(2,10)],'r');
plot([xa g(3)+Scale*Y(3,10)],[ya g(4)+Scale*Y(4,10)],'r');
plot([xb g(3)+Scale*Y(3,10)],[yb g(4)+Scale*Y(4,10)],'r');
plot([g(1)+Scale*Y(1,10) g(3)+Scale*Y(3,10)],[g(2)+Scale*Y(2,10) g(4)+Scale*Y(4,10)],'r');
plot([g(1)+Scale*Y(1,10) g(5)+Scale*Y(5,10)],[g(2)+Scale*Y(2,10) g(6)+Scale*Y(6,10)],'r');
plot([g(1)+Scale*Y(1,10) g(7)+Scale*Y(7,10)],[g(2)+Scale*Y(2,10) g(8)+Scale*Y(8,10)],'r');
plot([g(3)+Scale*Y(3,10) g(7)+Scale*Y(7,10)],[g(4)+Scale*Y(4,10) g(8)+Scale*Y(8,10)],'r');
plot([g(5)+Scale*Y(5,10) g(7)+Scale*Y(7,10)],[g(6)+Scale*Y(6,10) g(8)+Scale*Y(8,10)],'r');
plot([g(5)+Scale*Y(5,10) xg+Scale*Y(9,10)],[g(6)+Scale*Y(6,10) yg+Scale*Y(10,10)],'r');
plot([g(7)+Scale*Y(7,10) xg+Scale*Y(9,10)],[g(8)+Scale*Y(8,10) yg+Scale*Y(10,10)],'r');
axis equal; axis([-.5 3.5 -.5 1.5]);
title(sprintf('w = %f rad/s',sqrt(E(10,10))));




