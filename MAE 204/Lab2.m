%% Lab 2
%ethan foss

%% Measurements
% Distance from Inertial Axes:
d = 30;

% Link Measurements:
R1 = 23.5;
R2 = 12.5;
R23 = 24.5;
R3 = 10;
R34 = 21.5;
R4 = 10.5;
R5 = 8.9;
R6 = 9.5;

% Zero Position:
de1 = -46;
de2 = 7.5;
de3 = 14.6;

% Gripper Length:
lg = 13.4;

%% M Matrix:
M = [1 0 0 de1;0 1 0 de2-lg;0 0 1 de3;0 0 0 1];

%% Screw Axes:
w1 = [0 0 1]';
w2 = [0 -1 0]';
w3 = [0 -1 0]';
w4 = [0 -1 0]';
w5 = [0 0 -1]';
w6 = [0 -1 0]';

q1 = [0 d R1]';
q2 = q1+[0 -R2 0]';
q3 = q2+[-R23 0 0]';
q4 = q3+[-R34 R3-R4 0]';
q5 = q4+[0 0 -R5]';
q6 = q5+[0 -R6-lg 0]';

S1 = [w1;cross(w1,-q1)];
S2 = [w2;cross(w2,-q2)];
S3 = [w3;cross(w3,-q3)];
S4 = [w4;cross(w4,-q4)];
S5 = [w5;cross(w5,-q5)];
S6 = [w6;cross(w6,-q6)];

S = [S1 S2 S3 S4 S5 S6];

% Body Screw Axes:
for i = 1:6
    B(:,i) = Adjoint(M^-1)*S(:,i); % Body Screw Axes
end

%% Desired Positions:

pA = [-45;-6;7];
pB = [-45.6;-7.4;15];
pC = [-46;10;7];
pD = [-46;10;22];
pE = [-36;1;4.5];
pF = [-36;1;22];
pG = [-46;10;10.5];
% Desired Transformations:

TsdA = [M(1:3,1:3),pA;0 0 0 1];
TsdB = [M(1:3,1:3),pB;0 0 0 1];
TsdC = [M(1:3,1:3),pC;0 0 0 1];
TsdD = [M(1:3,1:3),pD;0 0 0 1];
TsdE = [M(1:3,1:3),pE;0 0 0 1];
TsdF = [M(1:3,1:3),pF;0 0 0 1];
TsdG = [M(1:3,1:3),pG;0 0 0 1];

letters = ["A"; "B"; "C"; "D"; "E"; "F"; "G"];
for i = 1:length(letters)
    disp("Desired End Effector Position of Point " + letters(i) + " in centimeters: " + mat2str(eval(["p" + letters(i)+"'"]),2) + " | Corresponding Transformation: " + mat2str(eval(["Tsd" + letters(i)]),2));
end

%% Inverse Kinematics:
ew = .001; %rad
ev = .1; %cm

thetaguess = [0;0;0;0;0;0];

[thetasA, successA,reportA] = IKinBodyIterations(B,M,TsdA,thetaguess,ew,ev);
[thetasB, successB,reportB] = IKinBodyIterations(B,M,TsdB,thetaguess,ew,ev);
[thetasC, successC,reportC] = IKinBodyIterations(B,M,TsdC,[0;0;45*pi/180;0;0;0],ew,ev);
[thetasD, successD,reportD] = IKinBodyIterations(B,M,TsdD,thetaguess,ew,ev);
[thetasE, successE,reportE] = IKinBodyIterations(B,M,TsdE,thetaguess,ew,ev);
[thetasF, successF,reportF] = IKinBodyIterations(B,M,TsdF,thetaguess,ew,ev);
[thetasG, successG,reportG] = IKinBodyIterations(B,M,TsdG,thetaguess,ew,ev);

if successA && successB && successC && successD && successE && successF && successG
    disp('All IK successfully found');
    letters = ["A"; "B"; "C"; "D"; "E"; "F"; "G"];
    for i = 1:length(letters)
        disp("Joint angles of Point " + letters(i) + " in degrees: " + mat2str(180/pi*eval(["thetas" + letters(i)+"'"]),3));
    end
else
    letters = ["A"; "B"; "C"; "D"; "E"; "F"; "G"];
    for i = 1:length(letters)
        if ~eval(["success" + letters(i)])
            disp("IK of " + letters(i) + " not found");
        end
    end
end
thetas = [thetasA';thetasB';thetasC';thetasD';thetasE';thetasF';thetasG']*180/pi;
writematrix(thetas,'C:\Users\ethan\OneDrive\Documents\MATLAB\MAE 204\Lab2Thetas.csv');
