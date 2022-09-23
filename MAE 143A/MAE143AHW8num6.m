clear all; clc;
%% number 6
syms n
x1(n) = piecewise(0<=n<=25,n*exp(-.3*n),26<=n<=49,0);
x2(n) = piecewise(0<=n<=25,n*exp(-.3*n),26<=n<=99,0);
n1 = 0:49;
n2 = 0:99;

X1 = fft(double(x1(n1)));
X2 = fft(double(x2(n2)));

figure(1); hold on;
subplot(2,1,1);
stem(n1,abs(X1)/50);
title('Magnitude of X1')
subplot(2,1,2);
stem(n1,angle(X1));
title('Angle of X1')

figure(2); hold on;
subplot(2,1,1);
stem(n2,abs(X2)/100);
title('Magnitude of X2')
subplot(2,1,2);
stem(n2,angle(X2));
title('Angle of X2')

x3 = [single(x1(n1)) single(x1(n1))];
X3 = fft(x3);

figure(3); hold on;
subplot(2,1,1);
stem(n2,abs(X3)/100);
title('Magnitude of X3')
subplot(2,1,2);
stem(n2,angle(X3));
title('Angle of X3')