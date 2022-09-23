%% test NextState
dt = .01;
Robot = struct('j',zeros(5,1),'w',zeros(4,1),'c',zeros(3,1));

jVels = [0,0,0,0,0]';
u = [10,10,10,10]';
u = [-10,10,-10,10]';
u = [-10,10,10,-10]';

for i = dt:dt:1
    Robot = NextState(Robot,jVels,u,dt,12.3,12.3)
end