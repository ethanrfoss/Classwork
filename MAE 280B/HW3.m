%% HW3

%% Problem 1

Kp = 6;
Kd = 6;
Ki = 3;
tau = .01;

G = tf([1],[1 0 0]);
Cd = tf([Kd 0],[tau 1]);
Cpi = tf([Kp Ki],[1 0]);

T1 = minreal(G*Cpi/(1+G*(Cpi+Cd)));
T2 = minreal(G*(Cpi+Cd)/(1+G*(Cpi+Cd)));

figure(1);
step(T1);
title('Non-ordinary Controller');
figure(2);
step(T2);
title('Ordinary Controller');

%% Problem 2

A = [0 1 0;0 0 0;1 0 0];
B = [0;1;0];
C = [1 0 0];

Q = diag([200 10 1]); R = 1;
[F,S,CLP] = lqr(A,B,Q,R);

K = tf(F(1)) + tf([F(2) 0],[tau 1]) + tf([F(3)],[1 0]);
[num,den] = ss2tf(A,B,C,0);
G = tf(num,den);

T = minreal(G*K/(1+G*K));

figure(3);
step(T);


%% Problem 4

m1 = 5.8*10^-6;
m2 = 7*10^-6;
m3 = 5*10^-6;
d1 = 1*10^-5;
d2 = 1*10^-6;
d3 = 1*10^-6;
k2 = 2.4*10^-3;
k3 = 2.4*10^-3;

A = [0 1 0 0 0 0;
     -k2/m1 (-d1-d2)/m1 k2/m1 d2/m1 0 0;
     0 0 0 1 0 0;
     k2/m2 d2/m2 (-k2-k3)/m2 (-d2-d3)/m2 k3/m2 d3/m2;
     0 0 0 0 0 1;
     0 0 k3/m3 d3/m3 -k3/m3 -d3/m3];
 
B = [0;1/m1;0;0;0;0];

C = [[1 0 0] zeros(1,3)];

R = 1;
Q = diag([20,.002,.001,.001,.15,.001]);

[K,S,CLP] = lqr(A,B,Q,R);

figure(4);
step(ss(A-B*K,[1;0;0;0;0;0],C,0));

disp('Controller Gains: ');

K


