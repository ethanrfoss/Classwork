
function HW9

a = .3;
x(:,1) = [-4.7;.2];
T = 10;
dt = 1/10000;
t = [0:dt:T];

for i = 1:length(t)
    x(:,i+1) = RK4(x(:,i),dt,@xdot1);
end
figure; hold on;
plot(x(1,:),x(2,:));
xr = [min(x(1,:))-1:.01:max(x(1,:))+1];
plot(xr,sin(xr).^2);
plot(xr,exp(a*xr));
plot(xr,2*exp(2*a*xr));
end

function xdot = xdot1(x)

eta = .001;
a = .3;
xdot = [-x(1)*x(2);-1/eta*(x(2)-(sin(x(1)))^2)*(x(2)-exp(a*x(1)))*(x(2)-2*exp(2*a*x(1)))];

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