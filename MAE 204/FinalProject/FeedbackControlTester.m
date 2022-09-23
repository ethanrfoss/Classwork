%% FeedbackControl Tester

Tsed = [0 0 1 .5 ;0 1 0 0;-1 0 0 .5;0 0 0 1];

Tsednext = [0 0 1 .6;0 1 0 0;-1 0 0 .3;0 0 0 1];

Tse = [.17 0 .985 .387;0 1 0 0;-.985 0 .17 .57;0 0 0 1];

Kp = zeros(6,6);
Ki = zeros(6,6);

dt = .01;

[Ve,Xerr] = FeedbackControl(Tse,Tsed,Tsednext,Kp,Ki,dt)

Robot = struct('j',zeros(5,1),'w',zeros(4,1),'c',zeros(3,1)); % 5 joint angles(j), 4 wheel angles(w), 3 chassis configuration paramaters(c)
Robot.j = [0;0;.2;-1.6;0]; % Initial Robot joint angles
Robot.w = [0;0;0;0]; % Inital Wheel Angles
Robot.c = [0;0;0]; % Initial Chassis Configuration
