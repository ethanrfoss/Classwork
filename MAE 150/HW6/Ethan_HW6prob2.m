clear all; clc;

%% 2a
A = pi*.1016^2;
E = 210*10^9;
rho = 7700;
L = 3.6576;
Ldiag = L*sqrt(2);
m = rho*A*L;
mDiag = rho*A*Ldiag;

k = @(b) [(cosd(b))^2 cosd(b)*sind(b) -(cosd(b))^2 -cosd(b)*sind(b) ; cosd(b)*sind(b) (sind(b))^2 -cosd(b)*sind(b) -(sind(b))^2; -(cosd(b))^2 -cosd(b)*sind(b) (cosd(b))^2 cosd(b)*sind(b); -cosd(b)*sind(b) -(sind(b))^2 cosd(b)*sind(b) (sind(b))^2];
k1to2 = A*E*k(0)/L;
k1to12 = A*E*k(45)/Ldiag;
k2to3 = A*E*k(0)/L;
k2to12 = A*E*k(90)/L;
k3to4 = A*E*k(0)/L;
k3to11 = A*E*k(90)/L;
k3to12 = A*E*k(135)/Ldiag;
k4to5 = A*E*k(0)/L;
k4to9 = A*E*k(45)/Ldiag;
k4to10 = A*E*k(90)/L;
k4to11 = A*E*k(135)/Ldiag;
k5to6 = A*E*k(0)/L;
k5to8 = A*E*k(45)/Ldiag;
k5to9 = A*E*k(90)/L;
k6to7 = A*E*k(0)/L;
k6to8 = A*E*k(90)/L;
k7to8 = A*E*k(135)/Ldiag;
k8to9 = A*E*k(180)/L;
k9to10 = A*E*k(180)/L;
k10to11 = A*E*k(180)/L;
k11to12 = A*E*k(180)/L;

kg = zeros(24,24);
kg([1,2,3,4],[1,2,3,4]) = kg([1,2,3,4],[1,2,3,4]) + k1to2;
kg([1,2,23,24],[1,2,23,24]) = kg([1,2,23,24],[1,2,23,24]) + k1to12;
kg([3,4,5,6],[3,4,5,6]) = kg([3,4,5,6],[3,4,5,6]) + k2to3;
kg([3,4,23,24],[3,4,23,24]) = kg([3,4,23,24],[3,4,23,24]) + k2to12;
kg([5,6,7,8],[5,6,7,8]) = kg([5,6,7,8],[5,6,7,8]) + k3to4;
kg([5,6,21,22],[5,6,21,22]) = kg([5,6,21,22],[5,6,21,22]) + k3to11;
kg([5,6,23,24],[5,6,23,24]) = kg([5,6,23,24],[5,6,23,24]) + k3to12;
kg([7,8,9,10],[7,8,9,10]) = kg([7,8,9,10],[7,8,9,10]) + k4to5;
kg([7,8,17,18],[7,8,17,18]) = kg([7,8,17,18],[7,8,17,18]) + k4to9;
kg([7,8,19,20],[7,8,19,20]) = kg([7,8,19,20],[7,8,19,20]) + k4to10;
kg([7,8,21,22],[7,8,21,22]) = kg([7,8,21,22],[7,8,21,22]) + k4to11;
kg([9,10,11,12],[9,10,11,12]) = kg([9,10,11,12],[9,10,11,12]) + k5to6;
kg([9,10,15,16],[9,10,15,16]) = kg([9,10,15,16],[9,10,15,16]) + k5to8;
kg([9,10,17,18],[9,10,17,18]) = kg([9,10,17,18],[9,10,17,18]) + k5to9;
kg([11,12,13,14],[11,12,13,14]) = kg([11,12,13,14],[11,12,13,14]) + k6to7;
kg([11,12,15,16],[11,12,15,16]) = kg([11,12,15,16],[11,12,15,16]) + k6to8;
kg([13,14,15,16],[13,14,15,16]) = kg([13,14,15,16],[13,14,15,16]) + k7to8;
kg([15,16,17,18],[15,16,17,18]) = kg([15,16,17,18],[15,16,17,18]) + k8to9;
kg([17,18,19,20],[17,18,19,20]) = kg([17,18,19,20],[17,18,19,20]) + k9to10;
kg([19,20,21,22],[19,20,21,22]) = kg([19,20,21,22],[19,20,21,22]) + k10to11;
kg([21,22,23,24],[21,22,23,24]) = kg([21,22,23,24],[21,22,23,24]) + k11to12;

M = @(l) rho*A*l*diag(.5*ones(1,4));
Mstr = M(L);
Mdiag = M(Ldiag);

Mg = zeros(24,24);
Mg([1,2,3,4],[1,2,3,4]) = Mg([1,2,3,4],[1,2,3,4]) + Mstr;
Mg([1,2,23,24],[1,2,23,24]) = Mg([1,2,23,24],[1,2,23,24]) + Mdiag;
Mg([3,4,5,6],[3,4,5,6]) = Mg([3,4,5,6],[3,4,5,6]) + Mstr;
Mg([3,4,23,24],[3,4,23,24]) = Mg([3,4,23,24],[3,4,23,24]) + Mstr;
Mg([5,6,7,8],[5,6,7,8]) = Mg([5,6,7,8],[5,6,7,8]) + Mstr;
Mg([5,6,21,22],[5,6,21,22]) = Mg([5,6,21,22],[5,6,21,22]) + Mstr;
Mg([5,6,23,24],[5,6,23,24]) = Mg([5,6,23,24],[5,6,23,24]) + Mdiag;
Mg([7,8,9,10],[7,8,9,10]) = Mg([7,8,9,10],[7,8,9,10]) + Mstr;
Mg([7,8,17,18],[7,8,17,18]) = Mg([7,8,17,18],[7,8,17,18]) + Mdiag;
Mg([7,8,19,20],[7,8,19,20]) = Mg([7,8,19,20],[7,8,19,20]) + Mstr;
Mg([7,8,21,22],[7,8,21,22]) = Mg([7,8,21,22],[7,8,21,22]) + Mdiag;
Mg([9,10,11,12],[9,10,11,12]) = Mg([9,10,11,12],[9,10,11,12]) + Mstr;
Mg([9,10,15,16],[9,10,15,16]) = Mg([9,10,15,16],[9,10,15,16]) + Mdiag;
Mg([9,10,17,18],[9,10,17,18]) = Mg([9,10,17,18],[9,10,17,18]) + Mstr;
Mg([11,12,13,14],[11,12,13,14]) = Mg([11,12,13,14],[11,12,13,14]) + Mstr;
Mg([11,12,15,16],[11,12,15,16]) = Mg([11,12,15,16],[11,12,15,16]) + Mstr;
Mg([13,14,15,16],[13,14,15,16]) = Mg([13,14,15,16],[13,14,15,16]) +Mdiag;
Mg([15,16,17,18],[15,16,17,18]) = Mg([15,16,17,18],[15,16,17,18]) + Mstr;
Mg([17,18,19,20],[17,18,19,20]) = Mg([17,18,19,20],[17,18,19,20]) + Mstr;
Mg([19,20,21,22],[19,20,21,22]) = Mg([19,20,21,22],[19,20,21,22]) + Mstr;
Mg([21,22,23,24],[21,22,23,24]) = Mg([21,22,23,24],[21,22,23,24]) + Mstr;

trace(kg)
trace(Mg)

%% 2b
kr = kg;
kr(:,14) = [];
kr(14,:) = [];
kr(:,2) = [];
kr(2,:) = [];
kr(:,1) = [];
kr(1,:) = [];
Mr = Mg;
Mr(:,14) = [];
Mr(14,:) = [];
Mr(:,2) = [];
Mr(2,:) = [];
Mr(:,1) = [];
Mr(1,:) = [];

[eigvec,eig] = eig(Mr\kr);
shapevecs = eigvec;
frequencies = diag(sqrt(eig))/(2*pi);
lowestfreqs = frequencies([15,16,21])
mode1 = eigvec(:,15);
mode2 = eigvec(:,16);
mode3 = eigvec(:,21);
mode1 = [0;0;mode1(1:11);0;mode1(12:end)];
mode2 = [0;0;mode2(1:11);0;mode2(12:end)];
mode3 = [0;0;mode3(1:11);0;mode3(12:end)];
%% 2d
tn = L*[0,0,1,0,2,0,3,0,4,0,5,0,6,0,5,1,4,1,3,1,2,1,1,1]';
figure(1); hold on;
subplot(3,1,1); hold on; title('Truss with Eigenmode deformations');
p1 = plot(tn([1,3]),tn([2,4]),'-k');
plot(tn([1,23]),tn([2,24]),'-k');
plot(tn([3,5]),tn([4,6]),'-k');
plot(tn([3,23]),tn([4,24]),'-k');
plot(tn([5,7]),tn([6,8]),'-k');
plot(tn([5,21]),tn([6,22]),'-k');
plot(tn([5,23]),tn([6,24]),'-k');
plot(tn([7,9]),tn([8,10]),'-k');
plot(tn([7,17]),tn([8,18]),'-k');
plot(tn([7,19]),tn([8,20]),'-k');
plot(tn([7,21]),tn([8,22]),'-k');
plot(tn([9,11]),tn([10,12]),'-k');
plot(tn([9,15]),tn([10,16]),'-k');
plot(tn([9,17]),tn([10,18]),'-k');
plot(tn([11,13]),tn([12,14]),'-k');
plot(tn([11,15]),tn([12,16]),'-k');
plot(tn([13,15]),tn([14,16]),'-k');
plot(tn([15,17]),tn([16,18]),'-k');
plot(tn([17,19]),tn([18,20]),'-k');
plot(tn([19,21]),tn([20,22]),'-k');
plot(tn([21,23]),tn([22,24]),'-k');

tn = tn+3*mode1;
p2 = plot(tn([1,3]),tn([2,4]),'--r');
plot(tn([1,23]),tn([2,24]),'--r');
plot(tn([3,5]),tn([4,6]),'--r');
plot(tn([3,23]),tn([4,24]),'--r');
plot(tn([5,7]),tn([6,8]),'--r');
plot(tn([5,21]),tn([6,22]),'--r');
plot(tn([5,23]),tn([6,24]),'--r');
plot(tn([7,9]),tn([8,10]),'--r');
plot(tn([7,17]),tn([8,18]),'--r');
plot(tn([7,19]),tn([8,20]),'--r');
plot(tn([7,21]),tn([8,22]),'--r');
plot(tn([9,11]),tn([10,12]),'--r');
plot(tn([9,15]),tn([10,16]),'--r');
plot(tn([9,17]),tn([10,18]),'--r');
plot(tn([11,13]),tn([12,14]),'--r');
plot(tn([11,15]),tn([12,16]),'--r');
plot(tn([13,15]),tn([14,16]),'--r');
plot(tn([15,17]),tn([16,18]),'--r');
plot(tn([17,19]),tn([18,20]),'--r');
plot(tn([19,21]),tn([20,22]),'--r');
plot(tn([21,23]),tn([22,24]),'--r');
legend([p1,p2],'Original Truss','Deformed Truss');
ylabel('y[m]');

tn = L*[0,0,1,0,2,0,3,0,4,0,5,0,6,0,5,1,4,1,3,1,2,1,1,1]';
subplot(3,1,2); hold on;
plot(tn([1,3]),tn([2,4]),'-k');
plot(tn([1,23]),tn([2,24]),'-k');
plot(tn([3,5]),tn([4,6]),'-k');
plot(tn([3,23]),tn([4,24]),'-k');
plot(tn([5,7]),tn([6,8]),'-k');
plot(tn([5,21]),tn([6,22]),'-k');
plot(tn([5,23]),tn([6,24]),'-k');
plot(tn([7,9]),tn([8,10]),'-k');
plot(tn([7,17]),tn([8,18]),'-k');
plot(tn([7,19]),tn([8,20]),'-k');
plot(tn([7,21]),tn([8,22]),'-k');
plot(tn([9,11]),tn([10,12]),'-k');
plot(tn([9,15]),tn([10,16]),'-k');
plot(tn([9,17]),tn([10,18]),'-k');
plot(tn([11,13]),tn([12,14]),'-k');
plot(tn([11,15]),tn([12,16]),'-k');
plot(tn([13,15]),tn([14,16]),'-k');
plot(tn([15,17]),tn([16,18]),'-k');
plot(tn([17,19]),tn([18,20]),'-k');
plot(tn([19,21]),tn([20,22]),'-k');
plot(tn([21,23]),tn([22,24]),'-k');

tn = tn+3*mode2;
plot(tn([1,3]),tn([2,4]),'--r');
plot(tn([1,23]),tn([2,24]),'--r');
plot(tn([3,5]),tn([4,6]),'--r');
plot(tn([3,23]),tn([4,24]),'--r');
plot(tn([5,7]),tn([6,8]),'--r');
plot(tn([5,21]),tn([6,22]),'--r');
plot(tn([5,23]),tn([6,24]),'--r');
plot(tn([7,9]),tn([8,10]),'--r');
plot(tn([7,17]),tn([8,18]),'--r');
plot(tn([7,19]),tn([8,20]),'--r');
plot(tn([7,21]),tn([8,22]),'--r');
plot(tn([9,11]),tn([10,12]),'--r');
plot(tn([9,15]),tn([10,16]),'--r');
plot(tn([9,17]),tn([10,18]),'--r');
plot(tn([11,13]),tn([12,14]),'--r');
plot(tn([11,15]),tn([12,16]),'--r');
plot(tn([13,15]),tn([14,16]),'--r');
plot(tn([15,17]),tn([16,18]),'--r');
plot(tn([17,19]),tn([18,20]),'--r');
plot(tn([19,21]),tn([20,22]),'--r');
plot(tn([21,23]),tn([22,24]),'--r');
ylabel('y[m]');

tn = L*[0,0,1,0,2,0,3,0,4,0,5,0,6,0,5,1,4,1,3,1,2,1,1,1]';
subplot(3,1,3); hold on;
plot(tn([1,3]),tn([2,4]),'-k');
plot(tn([1,23]),tn([2,24]),'-k');
plot(tn([3,5]),tn([4,6]),'-k');
plot(tn([3,23]),tn([4,24]),'-k');
plot(tn([5,7]),tn([6,8]),'-k');
plot(tn([5,21]),tn([6,22]),'-k');
plot(tn([5,23]),tn([6,24]),'-k');
plot(tn([7,9]),tn([8,10]),'-k');
plot(tn([7,17]),tn([8,18]),'-k');
plot(tn([7,19]),tn([8,20]),'-k');
plot(tn([7,21]),tn([8,22]),'-k');
plot(tn([9,11]),tn([10,12]),'-k');
plot(tn([9,15]),tn([10,16]),'-k');
plot(tn([9,17]),tn([10,18]),'-k');
plot(tn([11,13]),tn([12,14]),'-k');
plot(tn([11,15]),tn([12,16]),'-k');
plot(tn([13,15]),tn([14,16]),'-k');
plot(tn([15,17]),tn([16,18]),'-k');
plot(tn([17,19]),tn([18,20]),'-k');
plot(tn([19,21]),tn([20,22]),'-k');
plot(tn([21,23]),tn([22,24]),'-k');

tn = tn+3*mode3;
plot(tn([1,3]),tn([2,4]),'--r');
plot(tn([1,23]),tn([2,24]),'--r');
plot(tn([3,5]),tn([4,6]),'--r');
plot(tn([3,23]),tn([4,24]),'--r');
plot(tn([5,7]),tn([6,8]),'--r');
plot(tn([5,21]),tn([6,22]),'--r');
plot(tn([5,23]),tn([6,24]),'--r');
plot(tn([7,9]),tn([8,10]),'--r');
plot(tn([7,17]),tn([8,18]),'--r');
plot(tn([7,19]),tn([8,20]),'--r');
plot(tn([7,21]),tn([8,22]),'--r');
plot(tn([9,11]),tn([10,12]),'--r');
plot(tn([9,15]),tn([10,16]),'--r');
plot(tn([9,17]),tn([10,18]),'--r');
plot(tn([11,13]),tn([12,14]),'--r');
plot(tn([11,15]),tn([12,16]),'--r');
plot(tn([13,15]),tn([14,16]),'--r');
plot(tn([15,17]),tn([16,18]),'--r');
plot(tn([17,19]),tn([18,20]),'--r');
plot(tn([19,21]),tn([20,22]),'--r');
plot(tn([21,23]),tn([22,24]),'--r');
xlabel('x[m]');
ylabel('y[m]');