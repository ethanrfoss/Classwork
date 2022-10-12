%% Analysis of RK2 and RK4 for Lorenz and Rossler Systems
function RungeKuttaAnalysis

% Time Step and Terminal Time:
Tl = 10;
Tr = 100;
dt = 10.^[-1:-.1:-5];

% Lorenz Parameters:
LorenzParams.sigma = 4;
LorenzParams.b = 1;
LorenzParams.r = 48;

% Rossler Parameters:
RosslerParams.a = .2;
RosslerParams.b = .2;
RosslerParams.c = 5.7;

% Initial Conditions:
x0 = [-5;-5;0];

[~,xl] = RK4(x0,@lorenz,LorenzParams,10^-6,Tl);
[~,xr] = RK4(x0,@rossler,RosslerParams,10^-6,Tr);

% Simulate Lorenz and Rossler systems for multiple time steps:
for i = 1:length(dt)
    
    [~,x4l{i}] = RK4(x0,@lorenz,LorenzParams,dt(i),Tl);
    [~,x2l{i}] = RK2(x0,@lorenz,LorenzParams,dt(i),Tl,1/2);
    
    [~,x4r{i}] = RK4(x0,@rossler,RosslerParams,dt(i),Tr);
    [~,x2r{i}] = RK2(x0,@rossler,RosslerParams,dt(i),Tr,1/2);
    
    e4l(i) = norm(xl(:,end)-x4l{i}(:,end));
    e2l(i) = norm(xl(:,end)-x2l{i}(:,end));
    
    e4r(i) = norm(xr(:,end)-x4r{i}(:,end));
    e2r(i) = norm(xr(:,end)-x2r{i}(:,end));
end

% Plotting:
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure(1); 
subplot(2,2,1); hold on; title('Lorenz System','interpreter','latex');
plot3(xl(1,:),xl(2,:),xl(3,:));
subplot(2,2,2); hold on; title('Rossler System','interpreter','latex');
plot3(xr(1,:),xr(2,:),xr(3,:));
subplot(2,2,3); hold on;
loglog(-log10(dt),log10(e4l));
loglog(-log10(dt),log10(e2l));
legend('RK4','RK2');  xlabel('Log10 of Time Step','interpreter','latex'); ylabel('Log10 of final error','interpreter','latex');
subplot(2,2,4); hold on;
loglog(-log10(dt),log10(e4r));
loglog(-log10(dt),log10(e2r));
legend('RK4','RK2');  xlabel('Log10 of Time Step','interpreter','latex'); ylabel('Log10 of final error','interpreter','latex');

end

%% RK2 Method
% Inputs: 
% x - current state
% f - function handle to calculate derivative
% Params - Function Paramaters for f function
% dt - time step
% c - RK2 parameter
% Outputs:
% x - next state
function [t,x] = RK2(x0,f,Params,dt,T,c)

x(:,1) = x0;
t = 0:dt:T;

for i = 1:length(t)

    f1 = f(x(:,i),Params);
    f2 = f(x(:,i) + c*dt*f1,Params);

    x(:,i+1) = x(:,i) + dt*((1-1/(2*c))*f1+1/(2*c)*f2);

end

end

%% RK4 Method
% Inputs: 
% x - current state
% f - function handle to calculate derivative
% Params - Function Paramaters for f function
% dt - time step
% Outputs:
% x - next state
function [t,x] = RK4(x0,f,Params,dt,T)

x(:,1) = x0;
t = 0:dt:T;

for i = 1:length(t)

    f1 = f(x(:,i),Params);
    f2 = f(x(:,i) + dt*f1/2,Params);
    f3 = f(x(:,i) + dt*f2/2,Params);
    f4 = f(x(:,i) + dt*f3,Params);

    x(:,i+1) = x(:,i) + dt*(f1/6 + (f2+f3)/3 + f4/6);

end

end

%% Lorenz System Equations
% Inputs:
% x - current state
% LorenzParams - Lorenz Equation Parameters
% Outputs:
% f - State derivative
function f = lorenz(x,LorenzParams)

f = [LorenzParams.sigma*(x(2)-x(1)) ; -x(2)-x(1)*x(3) ; -LorenzParams.b*x(3)+x(1)*x(2)-LorenzParams.b*LorenzParams.r];

end

%% Rossler System Equations
% Inputs:
% x - current state
% RosslerParams - Rossler Equation Parameters
% Outputs:
% f - State derivative
function f = rossler(x,RosslerParams)

f = [ -x(2)-x(3) ; x(1) + RosslerParams.a*x(2) ; RosslerParams.b + x(3)*(x(1) -RosslerParams.c)];

end