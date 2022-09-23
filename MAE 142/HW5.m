%% Hw3

%% 2
[2.8 .25;-.78 -.38]^-1*[(39000/530)/(.5*.00028656*(1.8*968)^2);.025]

%% 3 

A = [-.15 -187;.012 -.4];
B= [6.2;-1.3];
C = [0 1];

[b,a] = ss2tf(A,B,C,0);

%% 4

syms YB

A = [YB/750 -747.5;.016 -.31];
Y = 0:10:3000;
figure(1); hold on;
for i = 1:length(Y)
    YB = Y(i);
    E = eig(double(subs(A)));
    plot(E,'o');
end

YB = 232.5+50;
E = eig(A);
solve(E(1))
solve(E(2))

fplot(E(1),[0 100000])

%% 5

wn = sqrt(7096/876);
damp = 1051/(876*wn)
