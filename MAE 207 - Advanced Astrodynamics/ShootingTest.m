function ShootingTest
m1 = 5.9722*10^24; %Earth Mass(kg)
m2 = 7.34767*10^22; %Moon Mass(kg)

mu = (m1/m2+1)^(-1); %Jacobi Mu Non-Dimensional Constant

p0 = [.5;0;0]; v0 =[0;1.15715;0];
x0 = [p0;v0];
P = 7.2;

Cd = 3;
C = @(x) 2*((1-mu)./sqrt((x(:,1)+mu).^2+x(:,2).^2+x(:,3).^2)+mu./sqrt((x(:,1)-1+mu).^2+x(:,2).^2+x(:,3).^2))+x(:,1).^2+x(:,2).^2-(x(:,4).^2+x(:,5).^2+x(:,6).^2);

xd = [0;0];
x0S = [1];
xfS = [2,4];

tf = P/2;
maxiter = 500;
tol = 10^(-8);
i = 0;

figure(1); hold on;

while 1
    
    tspan = linspace(0,tf,1000);
    xphi0 = [x0; reshape(eye(6,6),36,1)];
    [t,xphi] = ode45(@(t,xphi)XSTMEOM(t,xphi,mu),tspan,xphi0,odeset('AbsTol',1e-12,'RelTol',1e-9));
    
    xf = xphi(end,1:6)';
    phif = reshape(xphi(end,7:end),6,6);
    
    xfdot = JacobiEOM(t,xf,mu);
    xe = xd-xf(xfS);
    
    K = [phif(xfS,x0S) xfdot(xfS)];
    Delta = K'*inv(K*K')*xe;
    
    x0(x0S) = x0(x0S)+Delta(1:end-1);
    tf = tf + Delta(end);
    x0(5) = sqrt(2*((1-mu)/(x0(1)+mu)+mu/(x0(1)+mu-1))+x0(1)^2-Cd);
    
    plot3(xphi(:,1),xphi(:,2),xphi(:,3));
    
    if norm(xe)<tol || i>maxiter
        break;
    end
    
    i = i + 1;
    XT{i} = xphi(:,1:6);
    
end

if i>maxiter
    disp('Shoting Algorithm Did Not Converge Before Maximum Iteration Limit');
else
    disp(sprintf('Shoting Algorithm Converged in %d iterations',i));
end

P = tf*2;

end

function dxphi = XSTMEOM(t,xphi,mu)

x = xphi(1:6);
phi = reshape(xphi(7:end),6,6);

rho13 = ((x(1)+mu)^2+x(2)^2+x(3)^2)^(3/2);
rho23 = ((x(1)+mu-1)^2+x(2)^2+x(3)^2)^(3/2);
xdot = [x(4)
        x(5)
        x(6)
        -(1-mu)*(x(1)+mu)/rho13-mu*(x(1)+mu-1)/rho23+2*x(5)+x(1)     
        -((1-mu)/rho13+mu/rho23)*x(2)-2*x(4)+x(2)
        -((1-mu)/rho13+mu/rho23)*x(3)];

rho1 = (x(1)+mu).^2+x(2).^2+x(3).^2;
rho2 = (x(1)-1+mu).^2+x(2).^2+x(3).^2;

Uxx = [(mu - 1)/rho1^(3/2) - mu/rho2^(3/2) + (3*mu*(2*mu + 2*x(1) - 2)*(mu + x(1) - 1))/(2*rho2^(5/2)) - (3*(2*mu + 2*x(1))*(mu + x(1))*(mu - 1))/(2*rho1^(5/2)) + 1,(3*x(2)*mu*(mu + x(1) - 1))/rho2^(5/2) - (3*x(2)*(mu + x(1))*(mu - 1))/rho1^(5/2),(3*mu*x(3)*(mu + x(1) - 1))/rho2^(5/2) - (3*x(3)*(mu + x(1))*(mu - 1))/rho1^(5/2)
       -x(2)*((3*(2*mu + 2*x(1))*(mu - 1))/(2*rho1^(5/2)) - (3*mu*(2*mu + 2*x(1) - 2))/(2*rho2^(5/2))),(mu - 1)/rho1^(3/2) - mu/rho2^(3/2) + x(2)*((3*x(2)*mu)/rho2^(5/2) - (3*x(2)*(mu - 1))/rho1^(5/2)) + 1,x(2)*((3*mu*x(3))/rho2^(5/2) - (3*x(3)*(mu - 1))/rho1^(5/2))
       -x(3)*((3*(2*mu + 2*x(1))*(mu - 1))/(2*rho1^(5/2)) - (3*mu*(2*mu + 2*x(1) - 2))/(2*rho2^(5/2))),x(3)*((3*x(2)*mu)/rho2^(5/2) - (3*x(2)*(mu - 1))/rho1^(5/2)), x(3)*((3*mu*x(3))/rho2^(5/2) - (3*x(3)*(mu - 1))/rho1^(5/2)) + (mu - 1)/rho1^(3/2) - mu/rho2^(3/2)];
   
A = [zeros(3,3),eye(3,3)
     Uxx,[0 2 0;-2 0 0;0 0 0]];

phidot = A*phi;

dxphi = [xdot; reshape(phidot,36,1)];

end

function xdot = JacobiEOM(t,x,mu)

rho13 = ((x(1)+mu)^2+x(2)^2+x(3)^2)^(3/2);
rho23 = ((x(1)+mu-1)^2+x(2)^2+x(3)^2)^(3/2);
xdot = [x(4)
        x(5)
        x(6)
        -(1-mu)*(x(1)+mu)/rho13-mu*(x(1)+mu-1)/rho23+2*x(5)+x(1)     
        -((1-mu)/rho13+mu/rho23)*x(2)-2*x(4)+x(2)
        -((1-mu)/rho13+mu/rho23)*x(3)];

end
