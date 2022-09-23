function[uk,xk]=Example211(T,uk)
%STEP0: initialize simulation system & derived parameters
s.h=0.01;
s.N=T/s.h;
s.mc=10;
t=[0:s.N]*s.h;
s.m1=0.2;
s.L1=1;
s.ell1=s.L1/2;
s.I1=s.m1*s.ell1^2/12;
s.m2=0.1;
s.L2=0.5;
s.ell2=s.L2/2;
s.I2=s.m2*s.ell2^2/12;
alpha=0.1;
s.B=[0;0;0;1;0;0];
s.Q=1/10000*diag([2 20 20 0 0 0]);
s.R=0;
%s.QT=diag([4 380 120 1 480 1200]);
s.QT=diag([20 100 60 1 1000 2400]);
if nargin<2
    uk=zeros(s.N+1,1);
end
s.x0=[0;pi;pi;0;0;0];
xk(1:6,1)=s.x0;
res=0;
for k=0:1000
    k
    u=uk(1);
    x=s.x0;
    J=0.25*s.h*(x'*s.Q*x+u'*s.R*u);
    c=.5;
    %STEP 1: march/save state(from t=0−>T), compute cost
    for n=1:s.N
        u=uk(n);
        f1=RHS(x,u,s);
        f2=RHS(x+s.h*f1/2,u,s);
        f3=RHS(x+s.h*f2/2,u,s);
        f4=RHS(x+s.h*f3,u,s);
        x=x+s.h*(f1/6+(f2+f3)/3+f4/6);
        xk(1:6,n+1)=x;
        u=uk(n+1);
        xk(7:9,n)=f1(4:6);
        if n==s.N
            c=.25;
        end
        J=J+c*s.h*(x'*s.Q*x+u'*s.R*u);
    end
    f1=RHS(x,u,s);
    xk(7:9,s.N+1)=f1(4:6);
    E=ComputeE(x,s);
    J=J+0.5*(x'*E'*s.QT*E*x);
    r=s.QT*E*x;
    g(s.N+1,1)=s.B'*r+s.R*uk(s.N+1);
    %STEPS 2&3: march adjoint(from t=T−>0), compute gradient
    for n=s.N:-1:1
        xh=(xk(:,n+1)+xk(:,n))/2;
        f1=RHSa(r,xk(:,n+1),s);
        f2=RHSa(r-s.h*f1/2,xh,s);
        f3=RHSa(r-s.h*f2/2,xh,s);
        f4=RHSa(r-s.h*f3,xk(:,n),s);
        r=r-s.h*(f1/6+(f2+f3)/3+f4/6);
        g(n,1)=s.B'*r+s.R*uk(n);
    end
    res1=res;res=g'*g;
    %STEPS 4&5: update u and repeat
    if(mod(k,4)==0||alpha<1e-4)
        pk=-g;
    else
        pk=-g+pk*res/res1;
    end
    %conjugate gradient
    h = figure(1);clf;    
    subplot(2,1,1);
    plot(t,xk(1,:),'r-',t,xk(2,:),'b-',t,xk(3,:),'g-');
    legend('x','\theta1','\theta2')
    
    subplot(2,1,2);
    plot(t,uk,'r--');
    [AA,AB,AC,JA,JB,JC]=Bracket(@ComputeJEx211,0,alpha,J,uk,pk,s);%find triplet
    [alpha,J]=Brent(@ComputeJEx211,AA,AB,AC,JA,JB,JC,1e-5,uk,pk,s);%refine triplet
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 0
        imwrite(imind,cm,[cd '\SavedGifs\' 'MPC2' '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
    else
        imwrite(imind,cm,[cd '\SavedGifs\' 'MPC2' '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
    end

    uk=uk+alpha*pk; %update uk
    pause(0.01);
    if abs(alpha)<1e-12
        break;
    end
end
% s.mc=1;
% for n=1:s.N+1%Computeukcorrespondingtodifferents.mctogivesamexk
%     E=ComputeE(xk(1:6,n),s);N=ComputeN(xk(1:6,n),0,s);uk(n,1)=s.B'*(E*xk(4:9,n)-N);
% end
end%functionExample211
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=RHS(x,u,s)
E=ComputeE(x,s);
N=ComputeN(x,u,s);
R=E\N;
end%functionRHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=RHSa(r,x,s)
E=ComputeE(x,s);
A=ComputeA(x,s);
R=-E'\(A'*r+s.Q*x(1:6));
end%functionRHSa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E=ComputeE(x,s)
I=eye(3);
Z=zeros(3);
E=[I Z;Z [s.mc+s.m1+s.m2, -s.m1*s.ell1*cos(x(2)), -s.m2*s.ell2*cos(x(3));-s.m1*s.ell1*cos(x(2)), s.I1+s.m1*s.ell1^2, 0; -s.m2*s.ell2*cos(x(3)), 0, s.I2+s.m2*s.ell2^2]];
end%functionComputeE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=ComputeN(x,u,s)
N=[x(4);x(5);x(6);-s.m1*s.ell1*sin(x(2))*x(5)^2-s.m2*s.ell2*sin(x(3))*x(6)^2+u; s.m1*9.8*s.ell1*sin(x(2)); s.m2*9.8*s.ell2*sin(x(3))];
end%functionComputeN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=ComputeA(x,s)
g=9.8;
a42=s.m1*s.ell1*(x(8)*sin(x(2))+x(5)^2*cos(x(2)));
a45=2*s.m1*s.ell1*x(5)*sin(x(2));
a43=s.m2*s.ell2*(x(9)*sin(x(3))+x(6)^2*cos(x(3)));
a46=2*s.m2*s.ell2*x(6)*sin(x(3));
a52=s.m1*s.ell1*(g*cos(x(2))-x(7)*sin(x(2)));
a63=s.m2*s.ell2*(g*cos(x(3))-x(7)*sin(x(3)));
A=[zeros(3) eye(3);0 -a42 -a43 0 -a45 -a46;0 a52 0 0 0 0;0 0 a63 0 0 0];
end%functionComputeA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J=ComputeJEx211(utrial,s)
x=s.x0;
u=utrial(1);
J=0.25*s.h*(x'*s.Q*x+u'*s.R*u);
c=.5;
for n=1:s.N
    u=utrial(n);
    if n==s.N
        c=.25;
    end
    f1=RHS(x,u,s);
    f2=RHS(x+s.h*f1/2,u,s);
    f3=RHS(x+s.h*f2/2,u,s);
    f4=RHS(x+s.h*f3,u,s);
    x=x+s.h*(f1/6+(f2+f3)/3+f4/6);
    J=J+c*s.h*(x'*s.Q*x+u'*s.R*u);
end
E=ComputeE(x,s);
J=J+0.5*(x'*E'*s.QT*E*x);
end%functionComputeJEx211