clear all; clc;

A = zeros(1,100000);
B = zeros(1,100000);
C = zeros(1,100000);
n=1:100000;
for i=1:length(n)
    A(i) = (78-55)*rand +55;
    B(i) = 90 + 4*sqrt(2)*erfinv(2*rand-1);
    C(i) = .25*(sqrt(A(i))+sqrt(B(i)))^2;
end
figure(1); hold on;
subplot(3,1,1);
histogram(A, 100);
xlabel('Grade(%)');
ylabel('Frequency(#)');
title('Uniform distribution between 55 and 78');

subplot(3,1,2);
histogram(B, 100);
xlabel('Grade(%)');
ylabel('Frequency(#)');
title('Normal distribution with mean 90% and standard deviation 4%');

subplot(3,1,3);
histogram(C, 100);
xlabel('Grade(%)');
ylabel('Frequency(#)');
title('Combination of A and B');

over84counter = 0;
for i=1:length(n)
    if C(i)>=84
        over84counter = over84counter+1;
    end
end

over84chance = 100*over84counter/100000