function RobotRunner
%% Function Description:
% Wrapper Function for Robot Runner. This function implements NextState,
% TrajectoryGenerator, and FeedbackControl for the youbot to determine the
% series of joint movements that causes the robot to perform the block
% movement. The function saves the robot configuration at each time step to
% a csv file. It also plots the position/orientation error and joint
% angles, as a function of time, and animates the trajectory of the robot.
clear global;
%% User Inputs:
% Initial Robot Configuration:
Robot = struct('j',zeros(5,1),'w',zeros(4,1),'c',zeros(3,1)); % 5 joint angles(j), 4 wheel angles(w), 3 chassis configuration paramaters(c)
Robot.j = [0;.798;-1.122;-.745;0]; % Initial Robot joint angles
Robot.w = [0;0;0;0]; % Inital Wheel Angles
Robot.c = [0;0;0]; % Initial Chassis Configuration

% Reference Initial State:
Tseinitialreference = [0 0 1 0;0 1 0 0;-1 0 0 .5;0 0 0 1];

% k:
k = 1;

% Proportional Constants:
Kp = .5*eye(6,6);

% Integral Constants:
Ki = 1*eye(6,6);

% Cube Initial State:
Tscinitial = [1 0 0 1;0 1 0 0;0 0 1 .025;0 0 0 1]; % Standard Task
%Tscinitial = [0 -1 0 0;1 0 0 1;0 0 1 .025;0 0 0 1]; % New Task

% Cube Final State:
Tscfinal = [0 1 0 0;-1 0 0 -1;0 0 1 .025;0 0 0 1];

% Joint Angle Limits
JLims = struct('j1',[-pi/2 pi/2],'j2',[-1.553 1.17],'j3',[-1.45 .5],'j4',[-1.142 .5],'j5',[-pi pi]); %Joint Limits
%JLims = struct('j1',[-pi pi],'j2',[-pi pi],'j3',[-pi pi],'j4',[-pi pi],'j5',[-pi pi]); %No joint Limits

%% Robot Transformations, Jacobian, Chassis Odometry, and Other Important Parameters:
% Zero Position Transformation:
Moe = [1 0 0 .033;0 1 0 0;0 0 1 .6546;0 0 0 1];
% Space Frame to Chassis Body Transformation:
Tsb = @(Robot) [cos(Robot.c(1)) -sin(Robot.c(1)) 0 Robot.c(2); sin(Robot.c(1)) cos(Robot.c(1)) 0 Robot.c(3); 0 0 1 .0963;0 0 0 1];
% Chassis to Arm Transformation
Tbo = [1 0 0 .1662;0 1 0 0;0 0 1 .0026;0 0 0 1];
% Body Screw Axes:
B = [0 0 1 0 .033 0;0 -1 0 -.5076 0 0;0 -1 0 -.3526 0 0;0 -1 0 -.2176 0 0;0 0 1 0 0 0]';
% Arm Transformation as a function of Arm angles
Toe = @(Robot) Moe*MatrixExp6(VecTose3(B(:,1))*Robot.j(1))*MatrixExp6(VecTose3(B(:,2))*Robot.j(2))*MatrixExp6(VecTose3(B(:,3))*Robot.j(3))*MatrixExp6(VecTose3(B(:,4))*Robot.j(4))*MatrixExp6(VecTose3(B(:,5))*Robot.j(5));
% Space to End Effector Transfromation:
Tse = @(Robot) Tsb(Robot)*Tbo*Toe(Robot);
% F Matrix For Chassis:
r = .0475; l = .47/2; w = .3/2; m = 4;
F = r/4*[-1/(l+w) 1/(l+w) 1/(l+w) -1/(l+w);1 1 1 1;-1 1 -1 1];
% Joint/Wheel Speed Limits:
jVelMax = 12.3; wVelMax = 12.3;

%% Body Jacobian
% F6:
F6 = [zeros(2,m); F; zeros(1,m)];
% Full Jacobian([Jbase Jarm]):
Je = @(Robot) [Adjoint(Toe(Robot)^-1*Tbo^-1)*F6 , JacobianBody(B,Robot.j)];

%% Generate Trajectory:
% Create Graps and Standoff References:
% Tcestandoff:
Tcestandoff = [0 0 1 -.035;0 1 0 0;-1 0 0 .25;0 0 0 1];
% Tcegrasp:
Tcegrasp = [-.5 0 sqrt(3)/2 0;0 1 0 0;-sqrt(3)/2 0 -.5 0;0 0 0 1];
% Tseinitialactual from specified joint angles:
Tseinitialactual = Tse(Robot);
% Calculate Trajectory:
[T,gripState] = TrajectoryGenerator(Tseinitialreference,Tscinitial,Tscfinal,Tcegrasp,Tcestandoff,k);

%% Feedback Control and Marching Next States:
% dt:
dt = .01;
% Initialize Tactual:
Tactual = cell(1,length(T));
Tactual{1} = Tseinitialactual;
% Create ConfigK and XerrK cells to store every Kth iter
ConfigK  = cell(1,ceil(length(T)/k));
XerrK  = cell(1,ceil(length(T)/k));
gripStateK = zeros(1,ceil(length(T)/k));
for i = 1:length(T)-1
    % Find Error and Twist with FeedbackControl:
    [Ve,Xerr] = FeedbackControl(Tactual{i},T{i},T{i+1},Kp,Ki,dt);
    % Evaluate Next State:
    [Robot] = NextState(Robot,dt,jVelMax,wVelMax,Je,JLims,Ve);
    % Next Actual Configuration:
    Tactual{i+1} = Tse(Robot);
    if mod(i-1,k) == 0
        ConfigK{(i+k-1)/k} = Robot;
        XerrK{(i+k-1)/k} = Xerr;
        gripStateK((i+k-1)/k) = gripState(i);
    end        
end

%% Save csv file of states and EE error:
% csv file:
csvMat = zeros(length(ConfigK),13);
for n = 1:length(ConfigK)-1
    csvMat(n,:) = [ConfigK{n}.c(1), ConfigK{n}.c(2), ConfigK{n}.c(3), ConfigK{n}.j(1), ConfigK{n}.j(2), ConfigK{n}.j(3), ConfigK{n}.j(4),ConfigK{n}.j(5), ConfigK{n}.w(1), ConfigK{n}.w(2), ConfigK{n}.w(3), ConfigK{n}.w(4), gripStateK(n)];
end
writematrix(csvMat,'C:\Users\ethan\OneDrive\Documents\MATLAB\MAE 204\FinalProject\YouBotConfig.csv');
%writematrix(csvMat,[cd '\YouBotConfig.csv']);
save('C:\Users\ethan\OneDrive\Documents\MATLAB\MAE 204\FinalProject\Xerr.mat','XerrK');
%save([cd '\Xerr.mat'],'XerrK');

%% Plotting:
figure; 
%% Plot Xerr:
subplot(2,2,1); hold on;
for i = 1:length(XerrK)-1
    perr(i) = sqrt(sum(XerrK{i}(4:end).^2));
    oerr(i) = sqrt(sum(XerrK{i}(1:3).^2));
end
plot((.01:.01:30-.01),perr);
plot((0.01:.01:30-.01),oerr);
legend('Position Error','Orientation Error');
title('Error Plot'); xlabel('Time[s]'); ylabel('Error[m/s,rad/s]');
hold off;

%% Plot Joint Angles:
subplot(2,2,3); hold on;
for i = 1:length(XerrK)-1
    j1(i) = ConfigK{i}.j(1);
    j2(i) = ConfigK{i}.j(2);
    j3(i) = ConfigK{i}.j(3);
    j4(i) = ConfigK{i}.j(4);
    j5(i) = ConfigK{i}.j(5);
end
plot((.01:.01:30-.01),j1*180/pi);
plot((.01:.01:30-.01),j2*180/pi);
plot((.01:.01:30-.01),j3*180/pi);
plot((.01:.01:30-.01),j4*180/pi);
plot((.01:.01:30-.01),j5*180/pi);
legend('Joint 1','Joint 2','Joint 3','Joint 4','Joint 5');
title('Joint Angles'); xlabel('Time[s]'); ylabel('Angle[deg]');
hold off;

%% Plot Chassis Position:
subplot(2,2,2); hold on; axis equal;
for i = 1:length(XerrK)-1
    x(i) = ConfigK{i}.c(2);
    y(i) = ConfigK{i}.c(3);
end
plot(x,y);
%% Plot Wheel Angles:
subplot(2,2,4); hold on;
for i = 1:length(XerrK)-1
    w1(i) = ConfigK{i}.w(1);
    w2(i) = ConfigK{i}.w(2);
    w3(i) = ConfigK{i}.w(3);
    w4(i) = ConfigK{i}.w(4);
end
plot((.01:.01:30-.01),w1*180/pi);
plot((.01:.01:30-.01),w2*180/pi);
plot((.01:.01:30-.01),w3*180/pi);
plot((.01:.01:30-.01),w4*180/pi);
legend('Wheel 1','Wheel 2','Wheel 3','Wheel 4');
title('Wheel Angles'); xlabel('Time[s]'); ylabel('Angle[deg]');
hold off;
%% Plot Transforms
figure; hold on; axis equal; view(45,45);
xlabel('X[dm]');ylabel('Y[dm]');zlabel('Z[dm]');
anT = title('Transformation Animation(t = 0s)');
% Space Frame:
qs = quiver3(0,0,0,1,0,0,1,'k');
quiver3(0,0,0,0,1,0,1,'k');
quiver3(0,0,0,0,0,1,1,'k');
% Tseinital:
qsei = quiver3(Tseinitialreference(1,4)*10,Tseinitialreference(2,4)*10,Tseinitialreference(3,4)*10,Tseinitialreference(1,1),Tseinitialreference(2,1),Tseinitialreference(3,1),1,'g');
quiver3(Tseinitialreference(1,4)*10,Tseinitialreference(2,4)*10,Tseinitialreference(3,4)*10,Tseinitialreference(1,2),Tseinitialreference(2,2),Tseinitialreference(3,2),1,'g');
quiver3(Tseinitialreference(1,4)*10,Tseinitialreference(2,4)*10,Tseinitialreference(3,4)*10,Tseinitialreference(1,3),Tseinitialreference(2,3),Tseinitialreference(3,3),1,'g');
% Tscinitial:
qsci = quiver3(Tscinitial(1,4)*10,Tscinitial(2,4)*10,Tscinitial(3,4)*10,Tscinitial(1,1),Tscinitial(2,1),Tscinitial(3,1),1,'y');
quiver3(Tscinitial(1,4)*10,Tscinitial(2,4)*10,Tscinitial(3,4)*10,Tscinitial(1,2),Tscinitial(2,2),Tscinitial(3,2),1,'y');
quiver3(Tscinitial(1,4)*10,Tscinitial(2,4)*10,Tscinitial(3,4)*10,Tscinitial(1,3),Tscinitial(2,3),Tscinitial(3,3),1,'y');
% Tscfinal:
qscf = quiver3(Tscfinal(1,4)*10,Tscfinal(2,4)*10,Tscfinal(3,4)*10,Tscfinal(1,1),Tscfinal(2,1),Tscfinal(3,1),1,'r');
quiver3(Tscfinal(1,4)*10,Tscfinal(2,4)*10,Tscfinal(3,4)*10,Tscfinal(1,2),Tscfinal(2,2),Tscfinal(3,2),1,'r');
quiver3(Tscfinal(1,4)*10,Tscfinal(2,4)*10,Tscfinal(3,4)*10,Tscfinal(1,3),Tscfinal(2,3),Tscfinal(3,3),1,'r');
% Cubes:
cubeX = [.25 .25 .25 .25 -.25 -.25 -.25 -.25]; cubeY = [.25 .25 -.25 -.25 .25 .25 -.25 -.25]; cubeZ = [.25 -.25 .25 -.25 .25 -.25 .25 -.25];
cubehull = convhull(cubeX,cubeY,cubeZ);
tC = trisurf(cubehull,cubeX+Tscinitial(1,4)*10,cubeY+Tscinitial(2,4)*10,cubeZ+Tscinitial(3,4)*10,'FaceColor','g','EdgeColor','g');
% Iterative plot:
i = 1;
qT1 = quiver3(Tactual{i}(1,4)*10,Tactual{i}(2,4)*10,Tactual{i}(3,4)*10,Tactual{i}(1,1),Tactual{i}(2,1),Tactual{i}(3,1),1,'b');
qT2 = quiver3(Tactual{i}(1,4)*10,Tactual{i}(2,4)*10,Tactual{i}(3,4)*10,Tactual{i}(1,2),Tactual{i}(2,2),Tactual{i}(3,2),1,'b');
qT3 = quiver3(Tactual{i}(1,4)*10,Tactual{i}(2,4)*10,Tactual{i}(3,4)*10,Tactual{i}(1,3),Tactual{i}(2,3),Tactual{i}(3,3),1,'b');
for i = 1:k*5:length(T)
    set(qT1,'XData',Tactual{i}(1,4)*10,'YData',Tactual{i}(2,4)*10,'ZData',Tactual{i}(3,4)*10,'UData',Tactual{i}(1,1),'VData',Tactual{i}(2,1),'WData',Tactual{i}(3,1));
    set(qT2,'XData',Tactual{i}(1,4)*10,'YData',Tactual{i}(2,4)*10,'ZData',Tactual{i}(3,4)*10,'UData',Tactual{i}(1,2),'VData',Tactual{i}(2,2),'WData',Tactual{i}(3,2));
    set(qT3,'XData',Tactual{i}(1,4)*10,'YData',Tactual{i}(2,4)*10,'ZData',Tactual{i}(3,4)*10,'UData',Tactual{i}(1,3),'VData',Tactual{i}(2,3),'WData',Tactual{i}(3,3));
    set(anT,'String',sprintf('Transformation Animation(t = %.2fs)',i*dt));
    drawnow;
    if gripState(i) == 1
        Tsc = Tactual{i}*(Tcegrasp)^-1;
        Tsc(1:3,4) = Tsc(1:3,4)*10;
        C = Tsc*[cubeX;cubeY;cubeZ;ones(1,length(cubeX))];
        set(tC,'XData',C(1,:),'YData',C(2,:),'ZData',C(3,:));
    end
end


