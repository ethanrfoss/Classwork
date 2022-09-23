%% Midterm

%% 1
%% 2
%% 3
roots([170.6 305.7 435.1 0 23.1])
%% 6
syms s c
c = 1;
A = [0 c;-c 0];
(s*eye(2,2)-A)^-1
x0 = [1;0];
t = [0:.1:10];
for i = 1:length(t)
w(i) = ([1 0])*exp(A*t(i))*x0;
end
%% 8 

A = [-.16 -177.3;.019 -.36];
B = [10.72; -2.58];
C = [0 1];
D = 0;
[n,d] = ss2tf(A,B,C,D)
G = tf(n,d);
damp(G)

%% 9
Ga = @(B) -11*tf([1 .8],[1 B -.4]);
K = .25;
T = @(B) minreal(Ga(B)/(1-Ga(B)*K));
B = 0:.1:100;
figure; hold on;
for i = 1:length(B);
    plot(real(pole(T(B(i)))),imag(pole(T(B(i)))),'ro');
end

rlocus(tf([1 0],[1 2.75 1.8]))

%% 11
G = tf([1 .1],[1 -9 10 0]);
margin(G)

%% 12
Kq = 50;
N = 4/pi*1.4/(Kq*.02);
J = 6000;
G = tf([1],[J 0]);
C = 180/pi;
T = minreal(Kq*G*C/(1+Kq*G*C))
%% 13

G = tf([2],[1 10.5 5]);
T = @(K) minreal(G/(1+G*K));
Kt = 0:.01:10;

for i = 1:length(Kt)
    poles = pole(T(Kt(i)));
    if real(poles(1))<= -2 && real(poles(2))<= -2
        Kt 
        poles
        break;
    end
end