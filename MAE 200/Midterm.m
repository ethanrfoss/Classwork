%% Midterm MAE 200

%% 1
A = [5 4;4 3]; B = [-3 4;4 -5];

EA = diag([1 4 3]);
EB = diag([5 8 2]);

[~,S] =eig([1 3 1; 9 -4 8; 9 2 3]);

S*EA*S^-1

S*EB*S^-1

%% 2
A = [2 2 -2;5 1 -3;1 5 -3];

A = [2 2 2 2 -4;7 1 1 1 -5;1 7 1 1 -5;1 1 7 1 -5;1 1 1 7 -5];
A2 = A^2;
A3 = A^3;
A4 = A^4;
A5 = A^5;

%% 3

C = [4 1 2 3; 3 4 1 2;2 3 4 1;1 2 3 4];

D = [3 2 1 4;4 3 2 1;1 4 3 2;2 1 4 3];


%% 4
B = [1 0 0 0 0 0 0 0 0 0]';

Aa = [zeros(10,1) rand(10,9)];
ctrb(Aa,B)

Ab = rand(10,10);
Ab(2:10,1) = 0;
ctrb(Ab,B)

Ac = rand(10,10);
Ac(3:10,1) = 0;
C = ctrb(Ac,B)
rank(C)



