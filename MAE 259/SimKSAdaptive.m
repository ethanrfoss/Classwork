%% KS Simulation Adpative
% Simulates the Kuramoto Shivashinsky Equations with IMEX RK methods

function SimKSAdaptive

[Aim,bim,bimhat,Aex,bex,bexhat] = IEXRKCB3c; % Initialize Butcher Tableau

% Aim = [0 0 0 0;
%        4/15 4/15 0 0;
%        4/15 1/3 1/15 0;
%        4/15 1/3 7/30 1/6];
% 
% bim = [4/15;1/3;7/30;1/6];
% 
% Aex = [0 0 0 0;
%        8/15 0 0 0;
%        1/4 5/12 0 0;
%        1/4 0 3/4 0];
% 
% bex = [1/4;0;3/4;0];

% Allowable Error:
S.E = 10^-3;

% Sytem Parameters
S.L = 50; % Length
S.tf = 100; % Final Time
S.dt = .05; % Time Step
N = [128]; % For testing
%N = [2^11 2^13 2^15]; % REALLY Long Computation Time

% For Loop for Each Discretization Count
for i = 1:length(N)

    S.N = N(i);
    u0 = randn(N,1); % Initial Condition

    [t1,x1,u1,T1] = FullStorageIMEXRK(Aex,bex,bexhat,Aim,bim,bimhat,S,u0); % Full Storage Simulate
    disp(['N = ' num2str(N(i)) ' Full Storage Computation Time = '  num2str(T1)]);
    [t2,x2,u2,T2] = ThreeRegisterIMEXRK(Aex,bex,bexhat,Aim,bim,bimhat,S,u0); % Thre Register Storage Simulate
    disp(['N = ' num2str(N(i)) ' Three Register Computation Time = '  num2str(T2)]);
    [t3,x3,u3,T3] = TwoRegisterIMEXRK(Aex,bex,bexhat,Aim,bim,bimhat,S,u0); % Two Register Storage Simulate
    disp(['N = ' num2str(N(i)) ' Two Register Computation Time = '  num2str(T3)]);
    
    % Print Simulation Times
    save(sprintf('KSFull%d.mat',log2(N(i))),'t1','x1','u1','T1');
    save(sprintf('KSThree%d.mat',log2(N(i))),'t2','x2','u2','T2');
    save(sprintf('KSTwo%d.mat',log2(N(i))),'t3','x3','u3','T3');
    
    % Contour Plot
    figure; hold on;
    contour(x1,t1,u1',[.25 .75 1.25],'r-');  
    contour(x1,t1,u1',[-.25 -.75 -1.25],'b-.')

end

% Optional Animation of Results
figure; hold on; axis([0 S.L -1.5 1.5]);
p1 = plot(x1,u1(:,1));
p2 = plot(x2,u2(:,1));
p3 = plot(x3,u3(:,1));
tit = title(sprintf('Time = %1.1f [s]',t1(1)));
for i = 2:1:length(t1)

    set(p1,'XData',x1,'YData',u1(:,i));
    set(p2,'XData',x2,'YData',u2(:,i));
    set(p3,'XData',x3,'YData',u3(:,i));
    set(tit,'String',sprintf('Time = %1.1f [s]',t1(i)));

    drawnow;

end

end

%% KSg
% Function evaluation for nonlinear  non-stiff part of ODE
function g = KSg(uhat,N,kx)

uhat(fix(N/3)+1:end)=0; % dealias

r = NR_RFFTinv(uhat,N); % Take inverse

r = -r.^2; % square

g = 1i*kx.*NR_RFFT(r,N); % FFT

end

%% Full Storage IMEXRK
% Simulate KS with given butcher tableau, system parameters, and initial
% conditions
function [t,x,u,T] = FullStorageIMEXRK(Aex,bex,bexhat,Aim,bim,bimhat,S,u0)

% Integration Setup
s = length(Aex);
h(1) = S.dt;
t(1) = 0;
i = 1;

% Setup problem
dx = S.L/S.N;
x = (0:S.N-1)'*dx;
%u=zeros(S.N,length(t)); uhat = zeros(S.N/2,length(t));
u(:,1) = u0;
uhat = NR_RFFT(u(:,1),S.N);
uhat2 = uhat;
kx = (2*pi/S.L)*[0:S.N/2-1]';
Aop = kx.^2-kx.^4;

% Preallocate
f = zeros(S.N/2,s);
g = zeros(S.N/2,s);

% Timing
T = 0;

while true
    tic;
    % FLOP savings
    Aimdt = Aim*h(i);
    bimdt = bim*h(i);
    bimhatdt = bimhat*h(i);
    Aexdt = Aex*h(i);
    bexdt = bex*h(i);
    bexhatdt = bexhat*h(i);

    % Full Storage Algorithm Starts Here
    % k = 1
    y = uhat;
    uhat2 = uhat;
    f(:,1) = Aop.*y./(1-Aimdt(1,1)*Aop); % 4N FLOPs
    g(:,1) = KSg(y + Aimdt(1,1)*f(:,1),S.N,kx); % 2N FLOPs + 1 FFT

    for k = 2:s
        y = uhat + f(:,1:k-1)*Aimdt(k,1:k-1)' + g(:,1:k-1)*Aexdt(k,1:k-1)'; % 4Nk + 2N FLOPs s-1 times

        f(:,k) = Aop.*y./(1-Aimdt(k,k)*Aop); % 4N FLOPs s-1 times
        g(:,k) = KSg(y + Aimdt(k,k)*f(:,k),S.N,kx); % 2N Flops + 1 FFT s-1 times

    end
    
    uhat = uhat + f*bimdt + g*bexdt; % 4sN + 2N FLOPs
    uhat2 = uhat2 + f*bimhatdt + g*bexhatdt;

    % Full Storage Algorithm Ends Here 
    % Total Cost: 2Ns^2 + 14Ns - 4N FLOPs + s FFTs
    T = T + toc;

    u(:,i+1) = NR_RFFTinv(uhat,S.N); % Inv fft

    % Calculate Next Time Step and Break if Necessary:
    t(i+1) = t(i) + h(i);
    if t(i+1) >= S.tf
        break;
    end
    h(i+1) = min(h(i)*(S.E/norm(uhat2-uhat,inf))^(1/3),S.tf-t(i+1));

    i = i + 1;

end

figure; hold on;
title('Time Step of Full Register Scheme')
plot(1:length(h),h)

end

function [t,x,u,T] = TwoRegisterIMEXRK(Aex,bex,bexhat,Aim,bim,bimhat,S,u0)

% Integration Setup
t(1) = 0;
h(1) = S.dt;
s = length(Aex);
i = 1;

% Setup problem
dx = S.L/S.N;
x = (0:S.N-1)'*dx;
u(:,1) = u0;
uhat = NR_RFFT(u(:,1),S.N);
uhat2 = uhat;
kx = (2*pi/S.L)*[0:S.N/2-1]';
Aop = kx.^2-kx.^4;

% Timing
T = 0;

while true

    tic;
    % FLOP savings
    Aimdt = Aim*h(i);
    bimdt = bim*h(i);
    bimhatdt = bimhat*h(i);
    Aexdt = Aex*h(i);
    bexdt = bex*h(i);
    bexhatdt = bexhat*h(i);

    % Two Register Storage Algorithm Starts Here
    % k = 1
    uhat2 = uhat;
    y = uhat;

    y = y./(1-Aimdt(1,1)*Aop); % 3N FLOPs
    uhat = uhat + bimdt(1)*Aop.*y+bexdt(1)*KSg(y,S.N,kx); % 5N FLOPs + 1 FFT
    
    for k = 2:s

        y = uhat + (Aimdt(k,k-1)-bimdt(k-1))*Aop.*y + (Aexdt(k,k-1)-bexdt(k-1))*KSg(y,S.N,kx); % 7N FLOPs + 1 FFT s-1 times

        y = y./(1-Aimdt(k,k)*Aop); % 3N FLOPs s-1 times

        uhat = uhat + bimdt(k)*Aop.*y+bexdt(k)*KSg(y,S.N,kx); % 5N FLOPs + 1 FFT s-1 times
        uhat2 = uhat2 + bimhatdt(k)*Aop.*y+bexhatdt(k)*KSg(y,S.N,kx);


    end
    
    % Two Register Storage Algorithm Ends Here
    % Total Cost: 2Ns^2 + 14Ns - 4N FLOPs + s FFTs
    T = T + toc;

    u(:,i+1) = NR_RFFTinv(uhat,S.N);

    % Calculate Next Time Step and Break if Necessary:
    t(i+1) = t(i) + h(i);
    if t(i+1) >= S.tf
        break;
    end
    h(i+1) = min(h(i)*(S.E/norm(uhat2-uhat,inf))^(1/3),S.tf-t(i+1));

    i = i + 1;

end

figure; hold on;
title('Time Step of Two Register Scheme')
plot(1:length(h),h)

end

function [t,x,u,T] = ThreeRegisterIMEXRK(Aex,bex,bexhat,Aim,bim,bimhat,S,u0)

% Integration Setup
t(1) = 0;
h(1) = S.dt;
s = length(Aex);
i = 1;

% Setup problem
dx = S.L/S.N;
x = (0:S.N-1)'*dx;
u(:,1) = u0;
uhat = NR_RFFT(u(:,1),S.N);
uhat2 = uhat;
kx = (2*pi/S.L)*[0:S.N/2-1]';
Aop = kx.^2-kx.^4;

% Timing
T = 0;

while true
    tic;
    
    % FLOP savings
    Aimdt = Aim*h(i);
    bimdt = bim*h(i);
    bimhatdt = bimhat*h(i);
    Aexdt = Aex*h(i);
    bexdt = bex*h(i);
    bexhatdt = bexhat*h(i);

    % Three Register Storage Algorithm Starts Here

    % k = 1
    uhat2 = uhat;
    y = uhat;

    z = Aop.*y./(1-Aimdt(1,1)*Aop);
    y = KSg(y+Aimdt(1,1)*z,S.N,kx);
    uhat = uhat + bimdt(1)*z+bexdt(1)*y;

    for k = 2:s
        
        y = uhat + (Aimdt(k,k-1)-bimdt(k-1))*z + (Aexdt(k,k-1)-bexdt(k-1))*y;

        z = Aop.*y./(1-Aimdt(k,k)*Aop);
        y = KSg(y+Aimdt(k,k)*z,S.N,kx);

        uhat = uhat + bimdt(k)*z+bexdt(k)*y;

        uhat2 = uhat2 + bimhatdt(k)*z+bexhatdt(k)*y;


    end

    % Three Register Storage Algorithm Ednds Here

    T = T + toc;

    u(:,i+1) = NR_RFFTinv(uhat,S.N);

    % Calculate Next Time Step and Break if Necessary:
    t(i+1) = t(i) + h(i);
    if t(i+1) >= S.tf
        break;
    end
    h(i+1) = min(h(i)*(S.E/norm(uhat2-uhat,inf))^(1/3),S.tf-t(i+1));

    i = i + 1;

end

figure; hold on;
title('Time Step of Three Register Scheme')
plot(1:length(h),h)

end

%% FFT from NR
function [uh]=NR_RFFT(u,N)
% function [uh]=NR_RFFT(u,N)
% This routine was written by substituting RFFT2 into RFFT1 and simplifying.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 5.5.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap05">Chapter 5</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.
% See also RFFTinv, RFFT1, RFFT2.  Verify with: RFFTtest.

wh=NR_FFTdirect(u(1:2:N-1)+i*u(2:2:N),N/2,-1);
M=N/2+2;
for n=2:N/4
  uh(n,1)  =(wh(n)+conj(wh(M-n))-i*exp(-2*pi*i*(n-1)/N)*(wh(n)-conj(wh(M-n))))/4;
  uh(M-n,1)=(conj(wh(n))+wh(M-n)-i*exp( 2*pi*i*(n-1)/N)*(conj(wh(n))-wh(M-n)))/4;
end
uh(1,1)=(real(wh(1))+imag(wh(1)))/2 + i*(real(wh(1))-imag(wh(1)))/2;
uh(N/4+1,1)=(real(wh(N/4+1))-i*imag(wh(N/4+1)))/2;
end % function NR_RFFT

%% FFTinv from NR
function [u]=NR_RFFTinv(uh,N)
% function [u]=NR_RFFTinv(uh,N)
% This routine was written by inverting the steps of RFFT and doing them in reverse.
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 5.5.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap05">Chapter 5</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.
% See also NR_RFFT.  Verify with: NR_RFFTtest.

wh(1)=real(uh(1,1))+imag(uh(1,1)) + i*(real(uh(1,1))-imag(uh(1,1)));
wh(N/4+1)=(real(uh(N/4+1,1)) -i*imag(uh(N/4+1,1)))*2;
M=N/2+2;
for n=2:N/4
  wh(n)=uh(n,1)+conj(uh(M-n,1))+(uh(n,1)-conj(uh(M-n,1)))*i*exp(2*pi*i*(n-1)/N);
  wh(M-n)=conj(uh(n,1))+uh(M-n,1)+(conj(uh(n,1))-uh(M-n,1))*i*exp(-2*pi*i*(n-1)/N);
end
w=NR_FFTdirect(wh,N/2,1);
u(1:2:N-1,1)=real(w)';  u(2:2:N,1)=imag(w)';
end % function NR_RFFTinv

%% FFTdirect from NR
function x=NR_FFTdirect(x,N,g)
% function x=NR_FFTdirect(x,N,g)
% Compute the forward FFT (g=-1) or inverse FFT (g=1) of a vector x of order N=2^s.
% At each stage, defining Ns=2^stage, the elements divide into N/Ns groups, each with Ns
% elements.  Each group is split in half and combined as in Fig 5.2.
% The corresponding wavenumber vector is:  k=(2*pi/L)*[[0:N/2]';[-N/2+1:-1]'].
% See <a href="matlab:NRweb">Numerical Renaissance: simulation, optimization, & control</a>, Section 5.4.1.
% Part of <a href="matlab:help NRC">Numerical Renaissance Codebase 1.0</a>, <a href="matlab:help NRchap05">Chapter 5</a>; please read the <a href="matlab:help NRcopyleft">copyleft</a>.
% See also FFTrecursive, FFTnonreordered.  Verify with: FFTdirectTest.

s=log2(N);              % number of stages 
ind=bin2dec(fliplr(dec2bin([0:N-1]',s))); x=x(ind+1); % Put input into bit reversed order.
for stage=1:s           % For each stage...
   Ns=2^stage;          % Determine size of the transform computed at this stage.
   M=Ns/2; w=exp(g*2*pi*i/Ns);  % Calculate w
   for m=0:M-1          % Split each group in half to combine as in Fig 5.2:
      wm=w^m;           % Calculate the necessary power of w.
      for n=1:Ns:N      % Determine each pair of elements to be combined and combine them.
         a=m+n; b=M+a; x([a b])=x(a)+[1; -1]*x(b)*wm;
      end
   end
end
if g==-1, x=x/N; end              % Scale the forward version of the transform by 1/N.
end % function NR_FFTdirect

%% Butcehr Tableau Initialization
function [Aim,bim,bimhat,Aex,bex,bexhat] = IEXRKCB3c

aim21 = 0;
aim22 = 3375509829940/4525919076317;
aim32 = -11712383888607531889907/32694570495602105556248;
aim33 = 566138307881/912153721139;
b1 = 0;
b2 = 673488652607/2334033219546;
b3 = 493801219040/853653026979;
b4 = 184814777513/1389668723319;
aex21 = 3375509829940/4525919076317;
aex32 = 272778623835/1039454778728;
aex43 = 1660544566939/2334033219546;

bimhat1 = 0;
bimhat2 = 366319659506/1093160237145;
bimhat3 = 270096253287/480244073137;
bimhat4 = 104228367309/1017021570740;

bexhat1 = 449556814708/1155810555193;
bexhat2 = 0;
bexhat3 = 210901428686/1400818478499;
bexhat4 = 480175564215/1042748212601;

Aim = [0 0 0 0;
       aim21 aim22 0 0;
       b1 aim32 aim33 0;
       b1 b2 b3 b4];

bim = [b1;b2;b3;b4];
bimhat = [bimhat1;bimhat2;bimhat3;bimhat4];

Aex = [0 0 0 0;
       aex21 0 0 0;
       b1 aex32 0 0;
       b1 b2 aex43 0];

bex = [b1;b2;b3;b3];
bexhat = [bexhat1;bexhat2;bexhat3;bexhat4];

end
