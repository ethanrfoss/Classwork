function RobotNext = applyJointLimits(RobotNext,Robot,Je,Ve,dt)

j1lims = [-pi pi];
j2lims = [-1.177 1.553];
j3lims = [-1.45 -.15];
j4lims = [-1.142 -.15];
j5lims = [-pi pi];
recalculateAngles = false;

if RobotNext.j(1)<min(j1lims) || RobotNext.j(1)>max(j1lims)
    Je(:,5) = 0;  
    recalculateAngles = true;
end

if RobotNext.j(2)<min(j2lims) || RobotNext.j(2)>max(j2lims)
    Je(:,6) = 0;
    recalculateAngles = true;
end

if RobotNext.j(3)<min(j3lims) || RobotNext.j(3)>max(j3lims)
    Je(:,7) = 0;
    recalculateAngles = true;
end

if RobotNext.j(4)<min(j4lims) || RobotNext.j(4)>max(j4lims)
    Je(:,8) = 0;
    recalculateAngles = true;
end

if RobotNext.j(5)<min(j5lims) || RobotNext.j(5)>max(j5lims)
    Je(:,9) = 0;
    recalculateAngles = true;
end

if recalculateAngles
    Jinv = pinv(Je,.01);
    wVels = Jinv(1:4,:)*Ve; % Wheel Velocities
    jVels = Jinv(5:end,:)*Ve; % Joint Velocities
    RobotNext.w = Robot.w + wVels*dt;
    RobotNext.j = Robot.j + jVels*dt;
    r = .0475; l = .47/2; w = .3/2; m = 4;
    %H = [-l-w 1 -1; l+w 1 1;l+w 1 -1;-l-w 1 1];
    F = [-1/(l+w) 1/(l+w) 1/(l+w) -1/(l+w);1 1 1 1;-1 1 -1 1]*r/4;

    Vb = F*wVels;
    Robot.c = Robot.c + Vb*dt;
end

end