
function HW6


x(:,1) = [1;0;0];
T = 100;
dt = 1/1000;
t = [0:dt:T];

for i = 1:length(t)
    x(:,i+1) = RK4(x(:,i),dt,@xdot1);
end
figure
plot3(x(1,:),x(2,:),x(3,:));
end

function xdot = xdot1(x)

xdot = [x(2);-x(1)^3-x(2)^3+x(3)^3;-x(3)+x(2)];

end

function xdot = xdot2(x)

xdot = [x(2)+x(1)*x(3);-x(1)-x(2)+x(2)*x(3);-x(1)^2-x(2)^2];

end

function x = RK4(x,dt,xdot)

f1 = xdot(x);
f2 = xdot(x+f1*dt/2);
f3 = xdot(x+f2*dt/2);
f4 = xdot(x+f3*dt);

x = x + dt*(f1/6+(f2+f3)/3+f4/6);

end