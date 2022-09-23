% Initialize all constants.
clear all; mw=0.027; r=0.034; mb=0.180; l=0.0477; G=35.5; g=9.8;
Ib=2.63e-4; Im=3.6e-8; sbar=0.003; omega_f=1760; k=sbar/omega_f;
Iw=2*(mw*r^2/2 + G^2*Im);
c1=Iw+(mw+mb)*r^2; c2=mb*r*l; c3=2*G^2*k;
c4=2*G*sbar; c5=Ib+mb*l^2; c6=mb*g*l;
% Now set up equations of motion derived from mechanical and electrical
% properties using the state vector x=(phi theta phidot thetadot)^T
Ebar=[1 0 0 0; 0 1 0 0; 0 0 c1 c2; 0 0 c2 c5];
Abar=[0 0 1 0; 0 0 0 1; 0 0 -c3 c3; 0 c6 c3 -c3];
Bbar=[0; 0; c4; -c4]; A=Ebar\Abar; lambda=eig(A); B=Ebar\Bbar;
C=[0 0 1 -1; 0 1 0 0];
Con=[B A*B A*A*B A*A*A*B]; rank_Con=rank(Con)
Obs=[C; C*A; C*A*A; C*A*A*A]; rank_Obs=rank(Obs)
[U,Sig,V]=svd(Obs); format short; sigma_Obs=diag(Sig)

Ebar2 = Ebar(2:4,2:4);
Abar2 = Abar(2:4,2:4);
Bbar2 = Bbar(2:4);

A2=Ebar2\Abar2;
B2=Ebar2\Bbar2;
C2 = C(:,2:4);

Con2=[B2 A2*B2 A2*A2*B2];
Obs2=[C2; C2*A2; C2*A2*A2];

C3 = [0 1 -1;1 0 0;0 1 0];
Obs3=[C3; C3*A2; C3*A2*A2];

C4 = [1 0 0];
Obs4=[C4; C4*A2; C4*A2*A2];

%% 3
Rre = [B2 A2*B2 A2*A2*B2];

Are = Rre^-1*A2*Rre;
Bre = Rre^-1*B2;
Cre = C2*Rre;

a0 = -Are(1,end);
a1 = -Are(2,end);
a2 = -Are(3,end);

R1 = [1 0 0;a2 1 0;a1 a2 1];

Rc = R1';

Ac = Rc^-1*Are*Rc;
Bc = Rc^-1*Bre;
Cc = Cre*Rc;

deseig = -abs(eig(Ac));
k1 = a2+deseig(1)+deseig(2)+deseig(3);
k2 = a1 - deseig(1)*deseig(2)-deseig(1)*deseig(3)-deseig(2)*deseig(3);
k3 = a0+deseig(1)*deseig(2)*deseig(3);
K = [k1 k2 k3];
ABK = Ac+Bc*K;

%% 4
Rob=[C4; C4*A2; C4*A2*A2]^-1;

Aob = Rob^-1*A2*Rob;
Bob = Rob^-1*B2;
Cc = C4*Rob;


