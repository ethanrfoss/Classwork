clear all; clc;

%% 1b

Q = 30;
lnode = 1;
Eall = 200*10^9;

I1 = .1^4/12;
I2 = pi*.03^4/4;
I3 = pi*.02^4/4;

F1 = @(FR1) FR1;
M1 = @(MR1) MR1;
F2 = 0;
M2 = 0;
F3 = 0;
M3 = 0;
F4 = 0;
M4 = 0;

F5 = Q*lnode/2;
M5 = Q*(lnode)^2/12;
F6 = Q*lnode;
M6 = 0;
F7 = Q*lnode;
M7 = 0;
F8 = Q*lnode;
M8 = 0;
F9 = Q*lnode;
M9 = 0;
F10 = Q*lnode;
M10 = 0;
F11 = Q*lnode;
M11 = 0;
F12 = Q*lnode;
M12 = 0;

F13 = Q*lnode/2-100;
M13 = -Q*(lnode)^2/12;
F14 = 0;
M14 = 0;
F15 = 0;
M15 = 0;
F16 = 0;
M16 = 0;
F17 = 0;
M17 = -2000;

Fr = [F2,M2,F3,M3,F4,M4,F5,M5,F6,M6,F7,M7,F8,M8,F9,M9,F10,M10,F11,M11,F12,M12,F13,M13,F14,M14,F15,M15,F16,M16,F17,M17]';
    
K = @(E,I,L) [12*E*I/L^3 6*E*I/L^2 -12*E*I/L^3 6*E*I/L^2; 6*E*I/L^2 4*E*I/L -6*E*I/L^2 2*E*I/L; -12*E*I/L^3 -6*E*I/L^2 12*E*I/L^3 -6*E*I/L^2; 6*E*I/L^2 2*E*I/L -6*E*I/L^2 4*E*I/L];
K1to4 = K(Eall,I1,lnode);
K5to12 = K(Eall,I2,lnode);
K13to16 = K(Eall,I3,lnode);

Kg = zeros(34,34);
Kg([1:4],[1:4]) = Kg([1:4],[1:4]) + K1to4;
Kg([3:6],[3:6]) = Kg([3:6],[3:6]) + K1to4;
Kg([5:8],[5:8]) = Kg([5:8],[5:8]) + K1to4;
Kg([7:10],[7:10]) = Kg([7:10],[7:10]) + K1to4;

Kg([9:12],[9:12]) = Kg([9:12],[9:12]) + K5to12;
Kg([11:14],[11:14]) = Kg([11:14],[11:14]) + K5to12;
Kg([13:16],[13:16]) = Kg([13:16],[13:16]) + K5to12;
Kg([15:18],[15:18]) = Kg([15:18],[15:18]) + K5to12;
Kg([17:20],[17:20]) = Kg([17:20],[17:20]) + K5to12;
Kg([19:22],[19:22]) = Kg([19:22],[19:22]) + K5to12;
Kg([21:24],[21:24]) = Kg([21:24],[21:24]) + K5to12;
Kg([23:26],[23:26]) = Kg([23:26],[23:26]) + K5to12;

Kg([25:28],[25:28]) = Kg([25:28],[25:28]) + K13to16;
Kg([27:30],[27:30]) = Kg([27:30],[27:30]) + K13to16;
Kg([29:32],[29:32]) = Kg([29:32],[29:32]) + K13to16;
Kg([31:34],[31:34]) = Kg([31:34],[31:34]) + K13to16;

Kr = Kg;
Kr(:,2) = [];
Kr(2,:) = [];
Kr(:,1) = [];
Kr(1,:) = [];

norm(Kg,2)
cond(Kr)
%% 1c
ur = Kr\Fr;

%% 1d
ug = [0;0;ur];
Fe = Kg*ug;

%% 1e
de = ug([1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33]);
figure(1); hold on;
title('Length vs. Displacement');
xlabel('Length(m)');
ylabel('Displacement(m)');
plot(0:16,de);

%% 1f
thetae = ug([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34]);
figure(2); hold on;
title('Length vs. Angular Displacement');
xlabel('Length(m)');
ylabel('Angular Displacement(radians)');
plot(0:16,thetae);