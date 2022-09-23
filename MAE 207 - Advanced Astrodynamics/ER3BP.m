function ER3BP

%% System Parameters(Earth-Moon):
m1 = 5.9722*10^24; %Earth Mass(kg)
m2 = 7.34767*10^22; %Moon Mass(kg)

mu = (m1/m2+1)^(-1); %Jacobi Mu Non-Dimensional Constant

e = .0549;

% Lagrange Points of Earth Moon System:
L1 = [.83692,0];
L2 = [1.15568,0];
L3 = [-1.00506,0];
L4 = [.48785,.86603];
L5 = [.48785,-.86603];

%% Jacobi Synodic Frame Plot:
p0 = [.6;0;.2]; v0 =[0;.4;0]; f0 = 0;
x0 = [p0;v0;f0];
tspan = [0 100];

[t,x] = ode45(@(t,x)ElipticalEOM(t,x,e,mu),tspan,x0,odeset('AbsTol',1e-12,'RelTol',1e-9));

C = @(x) 2*((1-mu)./sqrt((x(:,1)+mu).^2+x(:,2).^2+x(:,3).^2)+mu./sqrt((x(:,1)-1+mu).^2+x(:,2).^2+x(:,3).^2))+x(:,1).^2+x(:,2).^2-(x(:,4).^2+x(:,5).^2+x(:,6).^2);
E = -C(x);

figure(1);
subplot(1,3,1); hold on; axis equal; grid on; axis([-1.4 1.4 -1.4 1.4]);
xlabel('$\xi$','interpreter','latex'); ylabel('$\eta$','interpreter','latex'); %%title('Circular-Planar Synodic Frame Three-Body Problem','interpreter','latex','fontsize',20);
plot([-1.4 1.4],[0 0],'--k');
plot([0 0],[-1.4 1.4],'--k');
plot(-mu,0,'k.','MarkerSize',40);
plot(1-mu,0,'k.','MarkerSize',20);
plot(0,0,'k.');

plot(L1(1),L1(2),'b*');
plot(L2(1),L2(2),'b*');
plot(L3(1),L3(2),'b*');
plot(L4(1),L4(2),'b*');
plot(L5(1),L5(2),'b*');

xp = plot3(x(:,1),x(:,2),x(:,3));

subplot(1,3,2); hold on;
xlabel('Time','interpreter','latex'); ylabel('Jacobi Integral','interpreter','latex');%%title('Jacobi Integral','interpreter','latex','fontsize',35);
xE = plot(t,E);

subplot(1,3,3); 
xlabel('Time'); ylabel('Moon True Anomaly[rad]');
plot(t,x(:,7));

end

function xdot = ElipticalEOM(t,x,e,mu)

rho13 = ((x(1)+mu)^2+x(2)^2+x(3)^2)^(3/2);
rho23 = ((x(1)+mu-1)^2+x(2)^2+x(3)^2)^(3/2);

xdot = [x(4)
        x(5)
        x(6)
        (-(1-mu)*(x(1)+mu)/rho13-mu*(x(1)+mu-1)/rho23+x(1))/(1+e*cos(x(7)))+2*x(5)
        (-((1-mu)/rho13+mu/rho23)*x(2)+x(2))/(1+e*cos(x(7)))-2*x(4)
        (-((1-mu)/rho13+mu/rho23)*x(3)-e*cos(x(7)))/(1+e*cos(x(7)))
        (1+e*cos(x(7)))^2/(1-e^2)^(3/2)];

end