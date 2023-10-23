
function FinalProb2

k1 = 0;
k2 = 2;

S.f = k2;

x0 = [0;100;100];
T = 10;

[t,x] = ode45(@(t,x) xdot(S,x,t),[0 T],x0,odeset('AbsTol',1e-12,'RelTol',1e-9));
figure;
plot3(x(:,1),x(:,2),x(:,3));
figure; hold on;
plot(t,sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2));
plot(t,(1+k2)^(1/2)*exp(-(k1)/(2*(1+k2))*t).*norm(x0));

end

function xdot = xdot(S,x,t)

xdot = [-S.f*x(1)+(1+S.f)*x(3);-x(2)+x(2)*x(3);-x(1)-x(2)-x(2)^2];

end