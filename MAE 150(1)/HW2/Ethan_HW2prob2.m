%% Uniform Distribution
% P - Load: 100+-3 N
% L - Length: .5+-.02 m
% E - Young's Modulus: 200 GPa
% I - Area Moment of Inertia: .01+-.01 m^4

P = zeros(1,5000);
L = zeros(1,5000);
E = 200*10^9;
I = zeros(1,5000);
deflection = zeros(1,5000);

for n=1:5000
    P(n) = 100 + (rand-.5)*6;
    L(n) = .5 + (rand-.5)*.04;
    I(n) = .1 + (rand-.5)*.02;
    deflection(n) = P(n)*(L(n))^3/(E*I(n));
end

figure(1);
subplot(2,2,1);
histogram(P,60);
title('Load Histogram');
xlabel('P(N)');
ylabel('Frequency(#)');

subplot(2,2,2);
histogram(L,60);
title('Length Histogram');
xlabel('L(m)');
ylabel('Frequency(#)');

subplot(2,2,3);
histogram(I,60);
title('Inertia Histogram');
xlabel('I(m^4)');
ylabel('Frequency(#)');

subplot(2,2,4);
histogram(deflection,60);
title('Deflection Histogram');
xlabel('Deflection(m)');
ylabel('Frequency(#)');

%% Normal Distribution(T=3*standard dev)
P = zeros(1,5000);
L = zeros(1,5000);
E = 200*10^9;
I = zeros(1,5000);
deflection = zeros(1,5000);

for n=1:5000
    P(n) = 100 + 1*sqrt(2)*erfinv(2*rand-1);
    L(n) = .5 + .02/3*sqrt(2)*erfinv(2*rand-1);
    I(n) = .1 + .01/3*sqrt(2)*erfinv(2*rand-1);
    deflection(n) = P(n)*(L(n))^3/(E*I(n));
end

figure(2);
subplot(2,2,1);
histogram(P,60);
title('Load Histogram');
xlabel('P(N)');
ylabel('Frequency(#)');

subplot(2,2,2);
histogram(L,60);
title('Length Histogram');
xlabel('L(m)');
ylabel('Frequency(#)');

subplot(2,2,3);
histogram(I,60);
title('Inertia Histogram');
xlabel('I(m^4)');
ylabel('Frequency(#)');

subplot(2,2,4);
histogram(deflection,60);
title('Deflection Histogram');
xlabel('Deflection(m)');
ylabel('Frequency(#)');