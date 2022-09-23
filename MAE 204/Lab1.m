
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

%% M Matrix:
M = [1 0 0 de1;0 1 0 de2;0 0 1 de3;0 0 0 1];

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
q6 = q5+[0 -R6 0]';

S1 = [w1;cross(w1,-q1)];
S2 = [w2;cross(w2,-q2)];
S3 = [w3;cross(w3,-q3)];
S4 = [w4;cross(w4,-q4)];
S5 = [w5;cross(w5,-q5)];
S6 = [w6;cross(w6,-q6)];

S = [S1 S2 S3 S4 S5 S6];

%% Joint Variables and Transformation Matrix:
theta = [-30,-45,75,15,0,0]*pi/180;

T = ArmTransformSpace(M,S(:,1:i),theta(1:i));

%% Plot Orientation
figure(1); hold on; axis equal;

q1 = quiver3(T(1,4),T(2,4),T(3,4),T(1,1),T(2,1),T(3,1),5,'k');
quiver3(T(1,4),T(2,4),T(3,4),T(1,2),T(2,2),T(3,2),5,'k');
quiver3(T(1,4),T(2,4),T(3,4),T(1,3),T(2,3),T(3,3),5,'k');

q2 = quiver3(0,0,0,1,0,0,5,'b');
quiver3(0,0,0,0,1,0,5,'b');
quiver3(0,0,0,0,0,1,5,'b');

p = legend([q1 q2],'End Effector Position and Orientation(scaled by 5)','Inertial Frame Axes(scaled by 5)');