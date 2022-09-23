function [Robot] = NextState(Robot,dt,jVelMax,wVelMax,Je,JLims,Ve)
%% Function Description:
% For this component, you will write a function called NextState that uses the kinematics of the youBot (see MR
% Exercise 13.33), your knowledge of velocity kinematics, and your knowledge of the Euler method to predict how the
% robot will move in a small timestep given its current configuration and velocity.
%% Inputs:
% Robot - Contains the Robot state variables
% jVels - Contains the joint angle velocities
% wVels - Contains the wheel angle velocities
% dt - time step
% jVelMax - max joint angle velocity
% wVelMax - max wheel angle velocity
%% Outputs:
% Robot - Contains the new Robot State Variables

%% Calculate Joint and Wheel Velocities:
J = Je(Robot);
%J(:,7) = 0; % Lock 3rd joint option
Jinv = pinv(J,.01);
wVels = Jinv(1:4,:)*Ve; % Wheel Velocities
jVels = Jinv(5:end,:)*Ve; % Joint Velocities

%% Apply Joint Velocity Limits:
for i = 1:length(wVels)
    if wVels(i)>wVelMax
        wVels(i) = wVelMax;
    elseif wVels(i)<-wVelMax
        wVels(i) = -wVelMax;
    end
end
for i = 1:length(jVels)
    if jVels(i)>jVelMax
        jVels(i) = jVelMax;
    elseif jVels(i)<-jVelMax
        jVels(i) = -jVelMax;
    end
end

%% Determine new joint and wheel angles:
RobotNext = Robot;
RobotNext.w = Robot.w + wVels*dt;
RobotNext.j = Robot.j + jVels*dt;

%% Check and Apply Joint Limits:
recalculateAngles = false;

if RobotNext.j(1)<min(JLims.j1) || RobotNext.j(1)>max(JLims.j1)
    J(:,5) = 0;  
    recalculateAngles = true;
end

if RobotNext.j(2)<min(JLims.j2) || RobotNext.j(2)>max(JLims.j2)
    J(:,6) = 0;
    recalculateAngles = true;
end

if RobotNext.j(3)<min(JLims.j3) || RobotNext.j(3)>max(JLims.j3)
    J(:,7) = 0;
    recalculateAngles = true;
end

if RobotNext.j(4)<min(JLims.j4) || RobotNext.j(4)>max(JLims.j4)
    J(:,8) = 0;
    recalculateAngles = true;
end

if RobotNext.j(5)<min(JLims.j5) || RobotNext.j(5)>max(JLims.j5)
    J(:,9) = 0;
    recalculateAngles = true;
end

if recalculateAngles
    Jinv = pinv(J,.01);
    wVels = Jinv(1:4,:)*Ve; % Wheel Velocities
    jVels = Jinv(5:end,:)*Ve; % Joint Velocities
    RobotNext.w = Robot.w + wVels*dt;
    RobotNext.j = Robot.j + jVels*dt;
end

Robot = RobotNext;

%% Determine chassis state through 4-Wheel Mechanum Kinematics:
r = .0475; l = .47/2; w = .3/2; m = 4;
%H = [-l-w 1 -1; l+w 1 1;l+w 1 -1;-l-w 1 1];
F = [-1/(l+w) 1/(l+w) 1/(l+w) -1/(l+w);1 1 1 1;-1 1 -1 1]*r/4;

Vb = F*wVels;
Robot.c = Robot.c + Vb*dt;

end