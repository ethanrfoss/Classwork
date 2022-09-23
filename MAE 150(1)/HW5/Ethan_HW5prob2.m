clear all;
clc;

%% 2a
P = 90*10^3;
E = 210*10^9;
A = 25/100^2;
k = @(b) [(cosd(b))^2 cosd(b)*sind(b) -(cosd(b))^2 -cosd(b)*sind(b) ; cosd(b)*sind(b) (sind(b))^2 -cosd(b)*sind(b) -(sind(b))^2; -(cosd(b))^2 -cosd(b)*sind(b) (cosd(b))^2 cosd(b)*sind(b); -cosd(b)*sind(b) -(sind(b))^2 cosd(b)*sind(b) (sind(b))^2];

k1 = A*E*k(45)/sqrt(2); %1-2
k2 = A*E*k(0)/1;%1-3
k3 = A*E*k(90)/1;%2-3
k4 = A*E*k(atand(.5/1))/sqrt(1^2+.5^2);%3-4
k5 = A*E*k(-atand(.5/1))/sqrt(1^2+.5^2);%2-4

ke = [ k1(1:2,1:2)+k2(1:2,1:2), k1(1:2,3:4), k2(1:2,3:4), zeros(2,2); k1(3:4,1:2), k1(3:4,3:4)+k3(1:2,1:2)+k5(1:2,1:2),k3(1:2,3:4),k5(1:2,3:4); k2(3:4,1:2), k3(3:4,1:2), k2(3:4,3:4)+k3(3:4,3:4)+k4(1:2,1:2), k4(1:2,3:4); zeros(2,2), k5(3:4,1:2), k4(3:4,1:2), k4(3:4,3:4)+k5(3:4,3:4)] %global matrix
%% 2b
Fp = [0 0 -P 0 0]';
kep = ke;
kep(:,8) = []; kep(:,7) = []; kep(:,2) = [];
kep(8,:) = []; kep(7,:) = []; kep(2,:) = [];
d = kep\Fp;
de = [d(1); 0;d(2);d(3);d(4);d(5);0;0] %displacement matrix
Fe = ke*de %Force matrix
R1 = Fe(2); % Reaction at 1
R4 = Fe(8); % Reaction at 4
%% 2c
p1 = [0,0];
p2 = [1,1];
p3 = [1,0];
p4 = [2,.5];
figure(1); hold on;
title('Original and Deformed Truss(Scaling: 200)');
xlabel('x-axis(m)');
ylabel('y-axis(m)');
plot1 = plot([p1(1),p2(1)],[p1(2),p2(2)],'-r'); %1
plot([p1(1),p3(1)],[p1(2),p3(2)],'-r'); %2
plot([p3(1),p2(1)],[p3(2),p2(2)],'-r'); %3
plot([p3(1),p4(1)],[p3(2),p4(2)],'-r'); %4
plot([p4(1),p2(1)],[p4(2),p2(2)],'-r'); %5
axis([-.5 2.5 -.5 1.5]);
d = d*200;
plot2 = plot([p1(1)+d(1),p2(1)+d(2)],[p1(2),p2(2)+d(3)],'--g'); %1
plot([p1(1)+d(1),p3(1)+d(4)],[p1(2),p3(2)+d(5)],'--g'); %2
plot([p3(1)+d(4),p2(1)+d(2)],[p3(2)+d(5),p2(2)+d(3)],'--g'); %3
plot([p3(1)+d(4),p4(1)],[p3(2)+d(5),p4(2)],'--g'); %4
plot([p4(1),p2(1)+d(2)],[p4(2),p2(2)+d(3)],'--g'); %5
legend([plot1,plot2],'Original Truss','Deformed Truss'); 