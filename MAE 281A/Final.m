
function Final

S.e = 4;
S.k = 3;

x0 = [5000;5000];
T = 1000;

[t,x] = ode45(@(t,x) xdot(S,x,t),[0 T],x0,odeset('AbsTol',1e-12,'RelTol',1e-9));

plot(x(:,1),x(:,2));
plot(t,sqrt(x(:,1).^2+x(:,2).^2));

end

function xdot = xdot(S,x,t)

xdot = [x(2);-x(1)-S.e*(S.k+sin(t))*x(2)];

end