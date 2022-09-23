%% 4.23
syms s
G = 1/2040/(s+240/2040);
K = 5;
vtarg = 3
figure(1); hold on;
xlabel('t');
ylabel('v1');
title("4.23");
fplot(ilaplace(vtarg/s*(G*K/(1+G*K))));

%% 4.24
W=2000;
figure(2);hold on;
title("4.24");
xlabel('t');
ylabel('v1');
fplot(ilaplace(vtarg/s*(G*K/(1+G*K))+G*W/s/(1+G*K)));
axis([0 500 -20 30])

%% 4.25

xtarg = 10;
figure(3);hold on;
title("4.25");
xlabel('t');
ylabel('x1');
fplot(ilaplace(xtarg/s*(G*K/s/(1+G*K/s))));
axis([0 500 -20 30])

%% 4.27
figure(4);hold on;
title("4.27");
xlabel('t');
ylabel('x1');
fplot(ilaplace(xtarg/s*(G*K/s/(s+G*K/s))+G*W/s/(s+G*K/s)));
axis([0 500 -20 30])

%% 4.28
syms ki kp;
K = 20 + 1/s;
figure(5);hold on;
xlabel('t');
ylabel('x1');
title("4.28");
fplot(ilaplace(xtarg/s*(G*K/(s+G*K))+G*W/s/(s+G*K)));
axis([0 500 -20 30])