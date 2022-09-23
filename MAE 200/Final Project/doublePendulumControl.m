function [t,xk,u,K] = doublePendulumControl(Q,R,s)
g = 9.8;

% X = icare(s.E^-1*s.A,s.E^-1*s.B,Q,R);
% K = -R^-1*s.B'*X;

X = icare(s.A,s.B,Q,R,0,s.E);
K = -R^-1*s.B'*X*s.E

T = 10; % simulate for 10 sec
s.h=0.01;
t=[0:s.h:T];
xk(1:6,1)=s.x0;

for n=1:length(t)
    x = xk(1:6,n);
    if t(n) > 80
        u(n) = 0;
    else
        u(n)=K*x;
    end
    f1=RHS(x,u(n),s);
    f2=RHS(x+s.h*f1/2,u(n),s);
    f3=RHS(x+s.h*f2/2,u(n),s);
    f4=RHS(x+s.h*f3,u(n),s);
    x=x+s.h*(f1/6+(f2+f3)/3+f4/6);
    xk(1:6,n+1)=x;
    xk(7:9,n)=f1(4:6);
    if x(2) > pi/2 || x(2) < -pi/2
        disp("Pendulum Fell after " + num2str(t(n)) + "sec, Stopping Simulation");
        tFall = t(n);
        t(n+2:end) = [];
        break;
    elseif x(3) > pi/2 || x(3) < -pi/2
        disp("Pendulum Fell after " + num2str(t(n)) + "sec, Stopping Simulation");
        tFall = t(n);
        t(n+2:end) = [];
        break;
    end
end
xk(:,end) = [];
end