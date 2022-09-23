%% Final Progect

%% Double Pendulum Model
s.mc=10;
s.m1=0.2;
s.L1=1;
s.ell1=s.L1/2;
s.I1=s.m1*s.ell1^2/12;
s.m2=0.1;
s.L2=0.5;
s.ell2=s.L2/2;
s.I2=s.m2*s.ell2^2/12;
I=eye(3);
Z=zeros(3);
s.E=ComputeE([0;0;0;0;0;0],s);
s.A = ComputeA([0;0;0;0;0;0;0;0;0],s);
s.B=[0;0;0;1;0;0];
s.h = .01;

%% 3a
[ukfinite,xkfinite] = Example211(3);
xkfinite(:,end)
%[uk,xk] = Example211(3,uk);

plotPend(xkfinite,s);

%% 3b
clear xk;
Q = diag([900 1000 500 0.1 60 10]);
%Q = diag([4000 40 10 .001 60 10]);
l = 10;
R = l^2;

s.x0=[0.1;0;0;0;0;0];
[t1,xk1,u1,~] = doublePendulumControl(Q,R,s);

s.x0=[0;0.1;0;0;0;0];
[t2,xk2,u2,~] = doublePendulumControl(Q,R,s);

s.x0=[0;0;0.1;0;0;0];
[t3,xk3,u3,~] = doublePendulumControl(Q,R,s);

s.x0 = xkfinite(1:6,end);
[tinfinite,xkinfinite,ukinfinite,s.K] = doublePendulumControl(Q,R,s);
t = [0:s.h:3 tinfinite];
xk = [xkfinite xkinfinite];
uk = [ukfinite' ukinfinite];

figure; 
subplot(6,3,1);
plot(t1,xk1(1,:));
title('x_{0} = (0.1,0,0,0,0,0)');
ylabel('x[m]','fontweight','bold');
subplot(6,3,4);
plot(t1,xk1(2,:));
ylabel('\theta_{1}[rad]','fontweight','bold');
subplot(6,3,7);
plot(t1,xk1(3,:));
ylabel('\theta_{2}[rad]','fontweight','bold');
subplot(6,3,10);
plot(t1,xk1(4,:));
ylabel('\partialx/\partialt[m/s]','fontweight','bold');
subplot(6,3,13);
plot(t1,xk1(5,:));
ylabel('\partial\theta_{1}/\partialt[rad/s]','fontweight','bold');
subplot(6,3,16);
plot(t1,xk1(6,:));
ylabel('\partial\theta_{2}/\partialt[rad/s]','fontweight','bold');
xlabel('Time[s]');

subplot(6,3,2);
plot(t2,xk2(1,:));
title('x_{0} = (0,0.1,0,0,0,0)');
subplot(6,3,5);
plot(t2,xk2(2,:));
subplot(6,3,8);
plot(t2,xk2(3,:));
subplot(6,3,11);
plot(t2,xk2(4,:));
subplot(6,3,14);
plot(t2,xk2(5,:));
subplot(6,3,17);
plot(t2,xk2(6,:));
xlabel('Time[s]');

subplot(6,3,3);
plot(t3,xk3(1,:));
title('x_{0} = (0,0,0.1,0,0,0)');
subplot(6,3,6);
plot(t3,xk3(2,:));
subplot(6,3,9);
plot(t3,xk3(3,:));
subplot(6,3,12);
plot(t3,xk3(4,:));
subplot(6,3,15);
plot(t3,xk3(5,:));
subplot(6,3,18);
plot(t3,xk3(6,:));
xlabel('Time[s]');

figure;
subplot(6,1,1);
plot(t3,xkinfinite(1,:));
title('x_{0} = (0,0,0.1,0,0,0)');
subplot(6,1,2);
plot(t3,xkinfinite(2,:));
subplot(6,1,3);
plot(t3,xkinfinite(3,:));
subplot(6,1,4);
plot(t3,xkinfinite(4,:));
subplot(6,1,5);
plot(t3,xkinfinite(5,:));
subplot(6,1,6);
plot(t3,xkinfinite(6,:));
xlabel('Time[s]');

eig(s.E^-1*s.A+s.E^-1*s.B*s.K)

plotPend(xk,s);
%% 3c
s.C = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0];
obsv(s.A*s.E^-1,s.C*s.E^-1);
Q = eye(6,6);
a = .1;
R = a^2*eye(3,3);

P = icare(s.A',s.C',Q,R,0,s.E);
s.L = -s.E'*P*s.C'*R^-1;

eig(s.E^-1*(s.A+s.L*s.C))


%% 3d
clear xk
clear xkerr
clear xerr
x0 = [.1;0;0;0;0;0];
xerr0 = [0.1;0;0;0;0;0];
xk(1:6,1) = x0;
xkerr(1:6,1) = xerr0;
T = 80; % simulate for 10 sec
s.h=0.01;
t=[0:s.h:T];
for n=1:length(t)
    x = xk(1:6,n);
    xerr = xkerr(1:6,n);
    if t(n) > 80
        u(n) = 0;
    else
        u(n)=s.K*(x-xerr);
    end
    f1=RHS(x,u(n),s);
    f2=RHS(x+s.h*f1/2,u(n),s);
    f3=RHS(x+s.h*f2/2,u(n),s);
    f4=RHS(x+s.h*f3,u(n),s);
    x=x+s.h*(f1/6+(f2+f3)/3+f4/6);
    xk(1:6,n+1)=x;
    xk(7:9,n)=f1(4:6);
    xkerr(1:6,n+1) = exp(s.E^-1*(s.A+s.L*s.C)*t(n))*xerr0;
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
xkerr(:,end) = [];

figure; hold on;
subplot(6,2,1);
plot(t1,xk(1,:));
title('x''');
ylabel('x[m]','fontweight','bold');
subplot(6,2,3);
plot(t1,xk(2,:));
ylabel('\theta_{1}[rad]','fontweight','bold');
subplot(6,2,5);
plot(t1,xk(3,:));
ylabel('\theta_{2}[rad]','fontweight','bold');
subplot(6,2,7);
plot(t1,xk(4,:));
ylabel('\partialx/\partialt[m/s]','fontweight','bold');
subplot(6,2,9);
plot(t1,xk(5,:));
ylabel('\partial\theta_{1}/\partialt[rad/s]','fontweight','bold');
subplot(6,2,11);
plot(t1,xk(6,:));
ylabel('\partial\theta_{2}/\partialt[rad/s]','fontweight','bold');
xlabel('Time[s]');

subplot(6,2,2);
plot(t1,xk(1,:)-xkerr(1,:));
title('$$\hat{x}$$','Interpreter','Latex','fontweight','bold');
ylabel('x[m]','fontweight','bold');
subplot(6,2,4);
plot(t1,xk(2,:)-xkerr(2,:));
ylabel('\theta_{1}[rad]','fontweight','bold');
subplot(6,2,6);
plot(t1,xk(3,:)-xkerr(3,:));
ylabel('\theta_{2}[rad]','fontweight','bold');
subplot(6,2,8);
plot(t1,xk(4,:)-xkerr(4,:));
ylabel('\partialx/\partialt[m/s]','fontweight','bold');
subplot(6,2,10);
plot(t1,xk(5,:)-xkerr(5,:));
ylabel('\partial\theta_{1}/\partialt[rad/s]','fontweight','bold');
subplot(6,2,12);
plot(t1,xk(6,:)-xkerr(6,:));
ylabel('\partial\theta_{2}/\partialt[rad/s]','fontweight','bold');
xlabel('Time[s]');
sgtitle('3d');
plotPend(xk,s);

%% 3e
clear Kk

Q = diag([900 1000 500 0.1 60 10]);
%Q = diag([4000 40 10 .001 60 10]);
l = 20;
R = l^2;
X = icare(ComputeA(xkfinite(:,1),s),s.B,Q,R,0,ComputeE(xkfinite(:,1),s));
%[X,K] = icare(ComputeA([0;0;0;0;0;0;0;0;0],s),s.B,Q,1,0,ComputeE([0;0;0;0;0],s));
s.h=0.01;
T = 3;
t=[0:s.h:T];
for n=1:length(t)
    E = ComputeE(xkfinite(1:6,n),s);
    f1=RHSX(xkfinite(:,n),s,X,Q,R);
    f2=RHSX(xkfinite(:,n),s,X+s.h*f1/2,Q,R);
    f3=RHSX(xkfinite(:,n),s,X+s.h*f2/2,Q,R);
    f4=RHSX(xkfinite(:,n),s,X+s.h*f3,Q,R);
    X=X+s.h*(f1/6+(f2+f3)/3+f4/6);
    Kk(n,1:6) = -R^-1*s.B'*X*E;    
end

figure; hold on;
title('3e');
plot(t,Kk(:,1:6))
xlabel('t[sec]'); ylabel('K'); hold off;
%% 3f
clear Lk

Q = eye(6,6);
a = 100;
R = a^2*eye(3,3);
P = icare(ComputeA(xkfinite(:,1),s)',s.C',Q,R,0,ComputeE(xkfinite(:,1),s));
s.h=0.01;
T = 3;
t=[0:s.h:T];
for n= 1:length(t)
    E = ComputeE(xkfinite(1:6,n),s);
    f1=RHSP(xkfinite(:,n),s,P,Q,R);
    f2=RHSP(xkfinite(:,n),s,P+s.h*f1/2,Q,R);
    f3=RHSP(xkfinite(:,n),s,P+s.h*f2/2,Q,R);
    f4=RHSP(xkfinite(:,n),s,P+s.h*f3,Q,R);
    P=P+s.h*(f1/6+(f2+f3)/3+f4/6);
    Lk{n} = -s.E'*P*s.C'*R^-1;
    Lk1(n,1:6,1:3) = -s.E'*P*s.C'*R^-1;
end

figure; hold on; title('3f');
plot(t,Lk1(:,1,1));
plot(t,Lk1(:,2,1));
plot(t,Lk1(:,3,1));
plot(t,Lk1(:,4,1));
plot(t,Lk1(:,5,1));
plot(t,Lk1(:,6,1));
plot(t,Lk1(:,1,2));
plot(t,Lk1(:,2,2));
plot(t,Lk1(:,3,2));
plot(t,Lk1(:,4,2));
plot(t,Lk1(:,5,2));
plot(t,Lk1(:,6,2));
plot(t,Lk1(:,1,3));
plot(t,Lk1(:,2,3));
plot(t,Lk1(:,3,3));
plot(t,Lk1(:,4,3));
plot(t,Lk1(:,5,3));
plot(t,Lk1(:,6,3));
xlabel('t[sec]'); ylabel('K'); hold off;

%% 3g
clear x xpert upert 
t= [0:s.h:80];
xpert(1:6,1) = [0.1;0;0;0;0;0];
x(1:6,1) = xkfinite(1:6,1)+xpert(1:6,1);
for n=1:length(t)
    if n>301
        Kk(n,:) = s.K;
        xbar = 0;
        ubar = 0;
        xpert(1:6,n) = x(1:6,n);
    else
        xpert(1:6,n) = x(1:6,n)-xkfinite(1:6,n);
        xbar = xkfinite(1:6,n);
        ubar = ukfinite(n);
    end
    upert(n) = Kk(n,:)*xpert(1:6,n);
    u = ubar+upert(n);
    f1=RHS(x(1:6,n),u,s);
    f2=RHS(x(1:6,n)+s.h*f1/2,u,s);
    f3=RHS(x(1:6,n)+s.h*f2/2,u,s);
    f4=RHS(x(1:6,n)+s.h*f3,u,s);
    x(1:6,n+1)=x(1:6,n)+s.h*(f1/6+(f2+f3)/3+f4/6);
    %xpert(1:6,n+1) = x(1:6,n+1)-xkfinite(1:6,n+1);
end
%plotPend(x,s)

figure; hold on;
subplot(6,1,1);
plot(t1,xpert(1,:));
title('3g');
ylabel('x[m]','fontweight','bold');
subplot(6,1,2);
plot(t1,xpert(2,:));
ylabel('\theta_{1}[rad]','fontweight','bold');
subplot(6,1,3);
plot(t1,xpert(3,:));
ylabel('\theta_{2}[rad]','fontweight','bold');
subplot(6,1,4);
plot(t1,xpert(4,:));
ylabel('\partialx/\partialt[m/s]','fontweight','bold');
subplot(6,1,5);
plot(t1,xpert(5,:));
ylabel('\partial\theta_{1}/\partialt[rad/s]','fontweight','bold');
subplot(6,1,6);
plot(t1,xpert(6,:));
ylabel('\partial\theta_{2}/\partialt[rad/s]','fontweight','bold');
xlabel('Time[s]');

%% 3h
clear xk
clear xkerr
clear xerr
xpert(1:6,1) = [.1;0;0;0;0;0];
xerr(1:6,1) = [0.1;0;0;0;0;0];
x(1:6,1) = xkfinite(1:6,1) + xpert(1:6,1);
T = 80; % simulate for 80 sec
s.h=0.01;
t=[0:s.h:T];
for n=1:length(t)
    if n>301
        Kk(n,:) = s.K;
        xbar = 0;
        ubar = 0;
        xpert(1:6,n) = x(1:6,n);
        Lk{n} = s.L;
    else
        xpert(1:6,n) = x(1:6,n)-xkfinite(1:6,n);
        xbar = xkfinite(1:6,n);
        ubar = ukfinite(n);
    end
    upert(n) = Kk(n,:)*(xpert(1:6,n)-xerr(1:6,n));
    u(n) = ubar + upert(n);
    f1=RHS(x(1:6,n),u(n),s);
    f2=RHS(x(1:6,n)+s.h*f1/2,u(n),s);
    f3=RHS(x(1:6,n)+s.h*f2/2,u(n),s);
    f4=RHS(x(1:6,n)+s.h*f3,u(n),s);
    x(1:6,n+1)=x(1:6,n)+s.h*(f1/6+(f2+f3)/3+f4/6);
    %Rk4 for error:
    f1=RHSerr(xerr(1:6,n),s,Lk{n});
    f2=RHSerr(xerr(1:6,n)+s.h*f1/2,s,Lk{n});
    f3=RHSerr(xerr(1:6,n)+s.h*f2/2,s,Lk{n});
    f4=RHSerr(xerr(1:6,n)+s.h*f3,s,Lk{n});
    xerr(1:6,n+1)=xerr(1:6,n)+s.h*(f1/6+(f2+f3)/3+f4/6);
end
x(:,end) = [];
xerr(:,end) = [];

figure; hold on;
subplot(6,2,1);
plot(t1,xpert(1,:));
title('x''');
ylabel('x[m]','fontweight','bold');
subplot(6,2,3);
plot(t1,xpert(2,:));
ylabel('\theta_{1}[rad]','fontweight','bold');
subplot(6,2,5);
plot(t1,xpert(3,:));
ylabel('\theta_{2}[rad]','fontweight','bold');
subplot(6,2,7);
plot(t1,xpert(4,:));
ylabel('\partialx/\partialt[m/s]','fontweight','bold');
subplot(6,2,9);
plot(t1,xpert(5,:));
ylabel('\partial\theta_{1}/\partialt[rad/s]','fontweight','bold');
subplot(6,2,11);
plot(t1,xpert(6,:));
ylabel('\partial\theta_{2}/\partialt[rad/s]','fontweight','bold');
xlabel('Time[s]');

subplot(6,2,2);
plot(t1,xpert(1,:)-xerr(1,:));
title('$$\hat{x}$$','Interpreter','Latex','fontweight','bold');
ylabel('x[m]','fontweight','bold');
subplot(6,2,4);
plot(t1,xpert(2,:)-xerr(2,:));
ylabel('\theta_{1}[rad]','fontweight','bold');
subplot(6,2,6);
plot(t1,xpert(3,:)-xerr(3,:));
ylabel('\theta_{2}[rad]','fontweight','bold');
subplot(6,2,8);
plot(t1,xpert(4,:)-xerr(4,:));
ylabel('\partialx/\partialt[m/s]','fontweight','bold');
subplot(6,2,10);
plot(t1,xpert(5,:)-xerr(5,:));
ylabel('\partial\theta_{1}/\partialt[rad/s]','fontweight','bold');
subplot(6,2,12);
plot(t1,xpert(6,:)-xerr(6,:));
ylabel('\partial\theta_{2}/\partialt[rad/s]','fontweight','bold');
xlabel('Time[s]');
sgtitle('3h');

plotPend(x,s,'FullSystem');

%% 5a
syms f(x) g(x) dVdx(x)
f(x) =  x-x^3;
g(x) = 1;
dVdx(x) = -2*(x^3-x-x*sqrt(x^4-2*x^2+2));
V(x) = int(dVdx)+1.1478;

dt = .01; T = 10; x(1) = -.5;
umax = .6; umin = -.6;
t = 0:dt:T;
for i = 1:length(t)
    Vx = V(x(i));
    dV = dVdx(x(i));
    Lfv = dV*f(x(i));
    Lgv = dV*g(x(i));
    sigma = sqrt(Lfv^2+x(i)^2*Lgv^2);
    A = double(Lgv); b = double(-Lfv-sigma); H = 1;
    u(i) = quadprog(H,double(f(x(i))),A,b); %,[],[],umin,umax);
    if u(i)>umax
        u(i) = umax;
    elseif u(i)<umin
        u(i) = umin;
    end
    dx = f(x(i)) + g(x(i))*u(i);
    xnew = double(x(i) + dx*dt);
    if abs(xnew-x(i))<10^-5
        break;
    end
    x(i+1) = xnew;
end
figure;
subplot(3,1,1);
plot(t(1:length(x)),x);
title('x vs t');
xlabel('t[sec]');
ylabel('x');
subplot(3,1,2);
plot(t(1:length(x)),u);
title('u vs t');
xlabel('t[sec]');
ylabel('u');
subplot(3,1,3);
plot(t(1:length(x)),V(x));
title('V vs t');
xlabel('t[sec]');
ylabel('V');
