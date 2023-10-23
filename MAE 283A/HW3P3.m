
%% Q1

N = 10000;
G0=tf([1 -.1 -.9],[1 -1.4 0.5],1); H0=tf([1 -.1],[1 -.9],1); F0 = tf([.72],[1 -.1 0],1);
u = randn(N,1); e = randn(N,1)*sqrt(.1);
y = lsim(G0,u) + lsim(H0,e);
yP = lsim((1-F0)*G0,u) + lsim(F0,y);

E = y-yP;
% tn = 10;
% RE = xcorr(E,E,tn,'biased');
% plot(-tn:tn,RE);
disp('Question 1: Variance of \eta: ');
var(E)


%% Q2
% 
% load('R333_STEP_data.mat');
% load('R333_FRF_data.mat');

N1 = 80; N2 = 80;
u = ones(length(y),1);
t = 0:Ts:(length(y)-1)*Ts;


for i = 2:N2+1
    R(:,i-1) = y(i:N1+i-1) - y(1:N1);
end

for i = 3:N2+2
    Rb(:,i-2) = y(i:N1+i-1) - y(2:N1+1);
end

[U,S,V] = svd(R);

s = diag(S);
figure(1); hold on; title('Singular Values of Step Realization'); ylabel('Magnitude'); xlabel('Number');
plot(1:20,s(1:20),'r*');

%n = input('What Rank?(Recomended 6): ');
n = 2;

R1 = U(:,1:n)*S(1:n,1:n)^.5;
R2 = S(1:n,1:n)^.5*V(:,1:n)';

R1L = S(1:n,1:n)^-.5*U(:,1:n)';
R2R = V(:,1:n)*S(1:n,1:n)^-.5;

D = y(1);
C = R1(1,:);
B = R2(:,1);
A = R1L*Rb*R2R;

ysim = lsim(ss(A,B,C,D,Ts),u,t);
figure(2); hold on; title('Simulated vs. Actual Step Response'); xlabel('Time[s]'); ylabel('Output');
plot(t,ysim);
plot(t,y); legend('Simulated','Actual');

[num den] = ss2tf(A,B,C,D,1);
num = num(num~=0);
theta0 = [num den(2:end)];
disp('theta0 from Step Respone Realization: ');
theta0

Ghat1 = tf(num,den,Ts);

%% Q2 Part 2

N = 513; 
w=linspace(0,pi/1,N)';
X=[exp(-[length(den)-length(num):length(den)-1]*1j.*w) -exp(-[1:length(den)-1]*1j.*w).*Gspa];
theta=[real(X);imag(X)]\[real(Gspa);imag(Gspa)];

% max_par_diff = 10;
% counter = 0;
% while max_par_diff>1e-8 && counter < 100
%     Weight=[ones(N,1) exp(-1j*w) exp(-2j*w) exp(-3j*w) exp(-4j*w) exp(-5j*w) exp(-6j*w)]*[1;theta(8:13)];
%     %Weight=abs(Weight).^4;
%     X=[Weight./ones(N,1) Weight./exp(-1j*w) Weight./exp(-2j*w) Weight./exp(-3j*w) Weight./exp(-4j*w) Weight./exp(-5j*w) Weight./exp(-6j*w) Weight./-Gspa.*exp(-1j*w) Weight./-Gspa.*exp(-2j*w) Weight./-Gspa.*exp(-3j*w) Weight./-Gspa.*exp(-4j*w) Weight./-Gspa.*exp(-5j*w) Weight./-Gspa.*exp(-6j*w)];
%     theta_new=[real(X);imag(X)]\[real(Weight.\Gspa);imag(Weight.\Gspa)];
%     counter=counter+1
%     max_par_diff=max(abs(theta-theta_new));
%     theta=theta_new;
% end

Ghat=tf(theta(1:length(num))',[1 theta(length(num)+1:end)'],1);
[mghat,pghat]=bode(Ghat,w);



figure(3);
subplot(2,1,1);
loglog(w,squeeze(mghat),'g',w,abs(Gspa),'r'); title('Bode plot of Curve Fit'); legend('Ghat','G_{spa}');
ylabel('Magnitude');
subplot(2,1,2);
semilogx(w,squeeze(pghat),'g',w,angle(Gspa)*180/pi,'r');
xlabel('Frequency[rad/s]'); ylabel('Phase');

disp('theta0 from Frequency Domain Curve Fit: ');
theta
