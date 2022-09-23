%% HW 5
% Ethan Foss, A16129635

%% 5.8:
syms L;
S = [0 1 0 0 0 0; 0 0 0 0 1 0; 0 0 1 2*L 0 0]'; % Screw axes

Js = jacobiansym(S)

%% 5.19:
S = [0 0 0 0 0 1;1 0 0 0 0 0;0 0 1 1 0 0;0 0 1 1/sqrt(2) 1/sqrt(2) 0;1/sqrt(2) -1/sqrt(2) 0 0 0 -1/sqrt(2);0 0 0 0 1 0]'; % Screw Axes

Js = jacobiansym(S);
Js1to4 = Js(:,1:4) % Part a

singular = rank(S) == length(S) % Part b

ci = S'*[0;1;-1;1;0;0] % Part Ci
cii = S'*[1;-1;0;1;0;-1] % Part Cii

%% 5.25:

W1 = 109; W2 = 82; L1 = 425; L2 = 392; H1 = 89; H2 = 95; % lengths
S = [0 0 1 0 0 0;0 1 0 -H1 0 0;0 1 0 -H1 0 L1;0 1 0 -H1 0 L1+L2;0 0 -1 -W1 L1+L2 0;0 1 0 H2-H1 0 L1+L2]'; % Screw Axes
J = jacobiansym(S);
theta1 = pi/2;theta2 = pi/2;theta3 = pi/2;theta4 = pi/2;theta5 = pi/2;theta6 = pi/2;
J = subs(J); % substitute thetas

Jw = J(1:3,:) % part a
Jv = J(4:6,:) % part a

Aw = Jw*Jw';
Av = Jv*Jv';

[Sw,Ew] = eig(Aw);
[Sv,Ev] = eig(Av);
Sv = double([Sv(:,1)/sqrt(sum(Sv(:,1).^2)),Sv(:,2)/sqrt(sum(Sv(:,2).^2)),Sv(:,3)/sqrt(sum(Sv(:,3).^2))]);

PrincipleAxesJw = real(double([sqrt(Ew) Sw]))
PrincipleAxesJv = real(double([sqrt(Ev) Sv]))

[St,Et] = eig(Aw^-1);
[Sf,Ef] = eig(Av^-1);
Sf = double([Sf(:,1)/sqrt(sum(Sf(:,1).^2)),Sf(:,2)/sqrt(sum(Sf(:,2).^2)),Sf(:,3)/sqrt(sum(Sf(:,3).^2))]);

PrincipleAxesJt = real(double([sqrt(Et) St]))
PrincipleAxesJf = real(double([sqrt(Ef) Sf]))

%% 6.9

B = [0 0 1 0 2 0;0 0 1 0 1 0]';
M = [1 0 0 2;0 1 0 0;0 0 1 0;0 0 0 1];
Tsd = [-.5 -.866 0 .366;.866 -.5 0 1.366; 0 0 1 0;0 0 0 1];
thetaguess = [0 30]'*pi/180;
ew = .001; ev = 10^-4;

[thetalist, success] = IKinBodyIterations(B, M, Tsd, thetaguess, ew, ev);


%% Part 2:

W1 = 109; W2 = 82; L1 = 425; L2 = 392; H1 = 89; H2 = 95; % Lengths(mm)
M = [-1 0 0 L1+L2;0 0 1 W1+W2; 0 1 0 H1-H2; 0 0 0 1]; % Zero Position M matrix
S = [0 0 1 0 0 0;0 1 0 -H1 0 0;0 1 0 -H1 0 L1;0 1 0 -H1 0 L1+L2;0 0 -1 -W1 L1+L2 0;0 1 0 H2-H1 0 L1+L2]'; % Screw axes

for i = 1:6
    B(:,i) = Adjoint(M^-1)*S(:,i); % Body Screw Axes
end

Tsd = [0 1 0 -500;0 0 -1 100;-1 0 0 100;0 0 0 1]; % Desired End effector position
ew = .001; %rad 
ev = .001*1000; % mm
thetaguess = [3,-1,2,-1,-.5,-1.5]'; % Theta Guess

[theta,success] = IKinBodyIterations(B, M, Tsd, thetaguess, ew, ev); % Determine thetas, create report

