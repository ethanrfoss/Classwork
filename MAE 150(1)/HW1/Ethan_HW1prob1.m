clear all; clc;

x = [0:.01:5];
yc = cosh(x);
y2 = 1+x.^2/factorial(2);
y4 = 1+x.^2/factorial(2)+x.^4/factorial(4)+ x.^6/factorial(6);
y10 = 1+x.^2/factorial(2)+x.^4/factorial(4)+ x.^6/factorial(6)+ x.^8/factorial(8)+ x.^10/factorial(10)+x.^12/factorial(12)+ x.^14/factorial(14)+x.^16/factorial(16)+ x.^18/factorial(18);
figure(1); hold on;
plot(x,yc);
plot(x,y2);
plot(x,y4);
plot(x,y10);
title("Problem 1");
xlabel("x-axis");
ylabel("y_axis");
legend("cosh(x)","taylor n=2","taylor n=4","taylor n=10");