
%% MAE 283A Rectilinear Positioning System Identification Final Project
% Name: Ethan Foss
% Email: erfoss@ucsd.edu
% Date: December 2022

%% Description of Project:
% Non-parametric and parametric analysis of given data set

%% System Description:
% Rectilinear mechanical System - An electromotor generates the torque or forces for the flexible mechanical system. 
% For identification purposes, the rectilinear system is being excited with a force that is generated by this electromotor 
% whereas the position of all masses are being measured. Position measurements are available for different configurations 
% of the flexible systems and depend on the number and size of masses. Both the input signal, given by the applied voltage 
% to the electromotor, and the output signals, given by position measurement encoder counts, are available for identification purposes.

%% Data Description: 
% input u  : voltage to electromotor (sine sweep from 1 Hz till 5Hz, over 37 seconds with amp. of 0.35 Volts)
% output y : linear position of first (left) cart
% sampling : 0.005304 sec.
% data     : 6964 points

%% Clear all, close all, clc:
clear all;
close all;
clc;

%% Load Data:
load('r_232_1.mat');
Ts = .005304;
N = 6964;
t = 0:Ts:Ts*N;

disp('System: 2-3-2 Configuration of Rectilinear Mechanical System.');
disp('Input: Voltage to Electromotor. Output: Linear Position of Left Cart.');
disp('Data: Sine Sweep Input from 1 to 5 Hz. 6964 Data Points at 0.005304s Sampling Time.');

%% Estimate and Plot Auto-correlation:
T = N/2;
Ru = xcorr(u,u,T,'biased');

figure(1); hold on; title('Auto-Correlation of Input Signal u'); 
xlabel('\tau'); ylabel('R_u(\tau)');
plot(-T:T,Ru);


disp('Auto-Correlation of Input Signal with no downsampling or other manipulation.');
pause;

%% Spectrum of Input signal:
Pu = 1/N*fft([Ru(T+1:2*T+1); Ru(1:T)]); %Compute fft of rearranged correlation to remove non-real values
w = linspace(0,pi/Ts,N/2+1); % Frequency from 0 to pi/Ts

figure(2); 
loglog(w,Pu(1:N/2+1)); % Ignore repeated part of Spectral Estimate
title('Spectrum of Input Signal'); 
xlabel('\omega[rad/s]'); ylabel('|\Phi_u(\omega)|'); 

disp('Spectrum of Input Signal From previous Auto-Correlation Signal');
pause;

%% Non-Parametric Estimate of G0:
Ryu = xcorr(y,u,T,'biased');
gamma = 1000; %Hanning window width
Ruw = Ru.*[zeros((2*T-gamma)/2,1); hamming(gamma+1); zeros((2*T-gamma)/2,1)]; % Apply hanning window to correlation
Ryuw = Ryu.*[zeros((2*T-gamma)/2,1); hamming(gamma+1); zeros((2*T-gamma)/2,1)];

Gf = fft([Ryuw(T+1:2*T+1); Ryuw(1:T)])./fft([Ruw(T+1:2*T+1); Ruw(1:T)]); % Compute fft of reorder wighted correlations

figure(3); sgtitle('Frequency Response of System');
subplot(2,1,1);
loglog(w,abs(Gf(1:N/2+1)));
ylabel('|G_0(\omega)|');
subplot(2,1,2);
semilogx(w,angle(Gf(1:N/2+1))*180/pi);
xlabel('\omega[rad/s]'); ylabel('\angle G_0(\omega)');

disp('Non-Parametric estimate of G0');
pause;

%% Downsample data by a factor of 10:
ud = u(1:10:end);
yd = y(1:10:end);
Nd = length(ud)-1;
Tsd = Ts*10;
td = t(1:10:end);

%% Retake correlation and Spectrums:
T = Nd/2;
Ru = xcorr(ud,ud,T,'biased');

figure(1); hold on;
plot(-T:T,Ru); legend('Normal','Downsampled');

Pu = 1/Nd*fft([Ru(T+1:2*T+1); Ru(1:T)]); 
w = linspace(0,pi/Tsd,Nd/2+1);

figure(2); hold on;
loglog(w,abs(Pu(1:Nd/2+1))); 
legend('Normal','Downsampled');

Ryu = xcorr(yd,ud,T,'biased');
gamma = 100; %Hanning window width
Ruw = Ru.*[zeros((2*T-gamma)/2,1); hamming(gamma+1); zeros((2*T-gamma)/2,1)];
Ryuw = Ryu.*[zeros((2*T-gamma)/2,1); hamming(gamma+1); zeros((2*T-gamma)/2,1)];

Gf = fft([Ryuw(T+1:2*T+1); Ryuw(1:T)])./fft([Ruw(T+1:2*T+1); Ruw(1:T)]);

figure(3); sgtitle('Frequency Response of System');
subplot(2,1,1); hold on;
loglog(w,abs(Gf(1:Nd/2+1)));
ylabel('|G_0(\omega)|');
legend('Normal','Downsampled');
subplot(2,1,2); hold on;
semilogx(w,angle(Gf(1:Nd/2+1))*180/pi);
xlabel('\omega[rad/s]'); ylabel('\angle G_0(\omega)');

disp('Downsampling of signals to remove unnecessary spectrum in high frequency. Recomputed Auto-Correlation and Spectrum.');
pause;

%% LS FIR Identification:

M = 60; % Estimate N/10 Coefficients

for i = M:-1:1
    PHI(:,M-i+1) = ud(i:Nd-M+i); % Construct PHI matrix
end
Y = yd(M:Nd);

thetaFIR = PHI\Y; % Linear Regression

ysim = filter(thetaFIR,1,ud);

% Calculate Covariance matrix of estimate
v = yd - ysim;
lambda = var(v);
PN = lambda*1/Nd*(PHI'*PHI)^-1;
SE = diag(sqrt(PN));

% Determine Confidence bounds from covariance
up = thetaFIR + SE*tinv(.99,Nd-M);
down = thetaFIR - SE*tinv(.99,Nd-M);
vup = yd - filter(up,1,ud);
vdown = yd - filter(down,1,ud);

figure(4); hold on;
title('FIR Coefficients with Confidence Bounds');
xlabel('k'); ylabel('g(k)');
plot(1:M,thetaFIR);
plot(1:M,down,'--r');
plot(1:M,up,'--r');
legend('FIR Coefficients','99% Coefficient Bounds');

disp('Now using Least Squares to determine FIR coefficients. Confidence Intervals are also plotted using the covariance matrix of the estimate.');
pause;

% Determine Confidence bounds on Reu(Ljung 16.65)
Reu = xcorr(v,ud,T,'unbiased');
Re = xcorr(v,v,T,'unbiased');
confBound = sqrt(sum(Re.*Ru)/Nd)*3;

figure(5); hold on;
title('Residual of FIR Model');
xlabel('\tau'); ylabel('R_{\epsilonu}(\tau)');
plot(0:2*M,Reu(T+1+(0:2*M)));
plot(0:2*M,-confBound*ones(1,length(0:2*M)),'--r');
plot(0:2*M,confBound*ones(1,length(0:2*M)),'--r');

disp('Correlation of Equation Error with the input. These residuals illustrate the poor quality of this FIR model.');
pause;


%% Realization Algorithm:

N1 = 100; N2 = 100;

% Construct Hankel Matrices for general I/O algorithm
for i = 2:N2+1
    HY(:,i-1) = yd(i:N1+i-1);
    HU(:,i-1) = ud(i:N1+i-1);
    
    HYb(:,i-1) = yd(i+1:N1+i);
    HUb(:,i-1) = ud(i+1:N1+i);
end

Q = HY*(eye(N2,N2)-HU'*pinv(HU*HU')*HU);
Qb = HYb*(eye(N2,N2)-HUb'*pinv(HUb*HUb')*HUb);
[U,S,V] = svd(Q);

figure(6); hold on; title('Singular Values of System');
plot(1:length(diag(S)),diag(S),'*r');

disp('Now performing a Generalized Realization Algorithm on the data. Singular values of system are plotted. A Third Order Model is chosen');
pause;

% Choose model order of 3 and calculate A,C
n = 3;
Q1 = U(:,1:n)*S(1:n,1:n)^.5;
Q2 = S(1:n,1:n)^.5*V(:,1:n)';
Q1L = S(1:n,1:n)^-.5*U(:,1:n)';
Q2R = V(:,1:n)*S(1:n,1:n)^-.5;

C = Q1(1,:);
A = Q1L*Qb*Q2R;

% Perform LS to determine B and D matrices
x(:,1) = zeros(n,1);
for i = 1:Nd
    x(:,i+1) = A'*x(:,i)+C'*ud(i)';
end
PHI = [x' ud];
Y = yd;
BD = PHI\Y;
B = BD(1:n);
D = BD(end);

[num den] = ss2tf(A,B,C,D,1);
Ghat = tf(num,den,1);

% Calculate impulse response coefs
gR(1) = D;
for i = 2:M
    gR(i) = C*A^(i-2)*B;
end

figure(7); hold on; title('Impulse Response Coefficients from Realization Algorithm'); 
xlabel('k'); ylabel('g(k)');
plot(1:M,gR);

disp('Impulse Response Coefficients are Plotted. Note the difference from the FIR model computed earlier.');
pause;

% Compute residuals and confidence bounds
ysim = lsim(Ghat,ud);
ess = yd-ysim;
Reu = xcorr(ess,ud,T,'unbiased');
Re = xcorr(ess,ess,T,'unbiased');
confBound = sqrt(sum(Re.*Ru)/Nd)*3;

figure(8); title('Residual Analysis of General Realization Algortihm');
xlabel('\tau'); ylabel('R_{\epsilonu}(\tau)'); 
hold on;
plot(0:2*n,Reu(T+1+(0:2*n)));
plot(0:2*n,confBound*ones(length(0:2*n)),'r--');
plot(0:2*n,-confBound*ones(length(0:2*n)),'r--');

disp('Residual of Realization Algorithm Model');
pause;

%% Prediction Error Model

% Order Determination for ARX model:
tDat = iddata(y,u,Ts);
NN = struc(1:10,1:10,1:10);
V = arxstruc(tDat(1:ceil(N/2)),tDat(ceil(N/2)+1:N),NN);
order = selstruc(V,0);

% Arx Model With best order:
Marx = arx(tDat,order);

% Frequency Domain estimation with same order:
fDat = idfrd(Gf(1:Nd/2+1),w,Ts);
Mfarx = arx(fDat,order);

% Box-Jenkins:
Mbj = bj(tDat,[3 2 2 3 1]);

% Plot residuals for each model
figure(9); 
subplot(3,1,1);
resid(Marx,tDat); title('Residual Analysis of ARX model');
subplot(3,1,2);
resid(Mfarx,tDat);  title('Residual Analysis of Frequency Fitting ARX Model');
subplot(3,1,3);
resid(Mbj,tDat);  title('Residual Analysis of Box-Jenkins Model');

disp('Residual plots for an ARX, Frequency Fitting ARX, and Box-Jenkins model plotted respectively. A third order system is determined to best fit the validation data.');
pause;

% Evaluate models with simulation and prediction
ysimarx = sim(Marx,u);
ysimbj = sim(Mbj,u);
figure(10);
subplot(2,1,1); hold on;
plot(t,y);
plot(t,ysimarx);
plot(t,ysimbj); legend('Actual','ARX Model','BJ Model');
title('Simulated Responses');
xlabel('Time[sec]');
subplot(2,1,2);
compare(tDat,Marx,'r',Mbj,'b',1);

disp('Simulated outputs and predicted outputs for ARX and Box Jenkins Model');




