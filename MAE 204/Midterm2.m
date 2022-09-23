%% Midterm 2, MAE 204

%% 1

Tsc = [0 1 0 .5;-1 0 0 1.6;0 0 1 .2;0 0 0 1];
Tcb = [0 -sqrt(3)/2 1/2 1;1 0 0 1;0 1/2 sqrt(3)/2 1;0 0 0 1];
Tsd = [.1667 -.6498 .7416 .7; .9832 .1667 -.0749 .3; -.0749 .7416 .6667 1.3;0 0 0 1];

Tbd = (Tsc*Tcb)^-1*Tsd
    
Tsb = Tsc*Tcb;

figure; hold on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

%S:
quiver3(0,0,0,1,0,0,1,'k');
quiver3(0,0,0,0,1,0,1,'k');
quiver3(0,0,0,0,0,1,1,'k');

%C:
quiver3(Tsc(1,4),Tsc(2,4),Tsc(3,4),Tsc(1,1),Tsc(2,1),Tsc(3,1),1,'g');
quiver3(Tsc(1,4),Tsc(2,4),Tsc(3,4),Tsc(1,2),Tsc(2,2),Tsc(3,2),1,'g');
quiver3(Tsc(1,4),Tsc(2,4),Tsc(3,4),Tsc(1,3),Tsc(2,3),Tsc(3,3),1,'g');

%B:
quiver3(Tsb(1,4),Tsb(2,4),Tsb(3,4),Tsb(1,1),Tsb(2,1),Tsb(3,1),1,'k');
quiver3(Tsb(1,4),Tsb(2,4),Tsb(3,4),Tsb(1,2),Tsb(2,2),Tsb(3,2),1,'k');
quiver3(Tsb(1,4),Tsb(2,4),Tsb(3,4),Tsb(1,3),Tsb(2,3),Tsb(3,3),1,'k');

%D:
quiver3(Tsd(1,4),Tsd(2,4),Tsd(3,4),Tsd(1,1),Tsd(2,1),Tsd(3,1),1,'b');
quiver3(Tsd(1,4),Tsd(2,4),Tsd(3,4),Tsd(1,2),Tsd(2,2),Tsd(3,2),1,'b');
quiver3(Tsd(1,4),Tsd(2,4),Tsd(3,4),Tsd(1,3),Tsd(2,3),Tsd(3,3),1,'b');

Vb = se3ToVec(MatrixLog6(Tbd));
Vs = Adjoint(Tsb)*Vb;
%Vb:
pd = Tbd(1:3,4);
ps = Tsb(1:3,1:3)*pd;
quiver3(Tsb(1,4),Tsb(2,4),Tsb(3,4),ps(1),ps(2),ps(3),1,'g');
quiver3(Tsb(1,4),Tsb(2,4),Tsb(3,4),Vb(4),Vb(5),Vb(6),1,'Color',[.5 0 .5]);
quiver3(Tsb(1,4),Tsb(2,4),Tsb(3,4),Vb(1),Vb(2),Vb(3),1,'r');

Vd = Adjoint(Tbd^-1)*Vb;
quiver3(Tsd(1,4),Tsd(2,4),Tsd(3,4),Vb(4),Vb(5),Vb(6),1,'Color',[.5 0 .5]);

Vd = Adjoint(Tbd^-1)*Vb;


Rbd = Tbd(1:3,1:3);
theta = acos(.5*(trace(Rbd)-1));

wh = 1/(2*sin(theta))*(Rbd-Rbd^-1);
w =[-wh(2,3);wh(1,3);-wh(1,2)]
%% 2:

syms L
L = 1;
B = [0 1 0 -L 0 -3*L;-1 0 0 0 -L 0;0 0 0 0 0 -1;0 0 -1 0 -3*L 0;0 1 0 L 0 -2*L;0 0 -1 0 -L 0]';
Jb = jacobiansymB(B);

M = [0 0 -1 2*L;-1 0 0 -3*L; 0 1 0 -L;0 0 0 1];
for i = 1:6
    S(:,i) = real(Adjoint(M)*B(:,i)); % Body Screw Axes
end
Js = jacobiansym(S);

theta1 = 0; theta2 = 0; theta3 = 0; theta4 = 0; theta5 = pi; theta6 = 0;
Jbs = subs(Jb)
Jss = subs(Js)
Jbs*[1;0;1;0;0;1]

Ftip = [0 0 -1 0 0 -1]';
Felb = [0 0 2 0 0 -2]';
M = [0 0 1 4*L;1 0 0 L;0 1 0 -L;0 0 0 1];
Ftips = [M(1:3,1:3) M(1:3,1:3);M(1:3,1:3) M(1:3,1:3)]*Ftip;
Felbs = Adjoint(M)*Felb;
t = Jbs'*Ftip+[Jbs(:,1:3)';zeros(3,6)]*Felb;
ts = Jss'*Ftips+[Jss(:,1:3)';zeros(3,6)]*Felbs;

for i = 1:6
    S(:,i) = real(Adjoint(M)*Jbs(:,i)); % Body Screw Axes
end

Moc = [1 0 0 2*L;0 1 0 -L; 0 0 1 -L;0 0 0 1];
Jm = jacobiansym(
%% 3:

S= [0 0 1 0 0 0;0 0 1 2 0 0;0 0 1 3 0 0;0 0 0 0 0 1]';
M = [1 0 0 0 ;0 1 0 3;0 0 1 0;0 0 0 1];
thetaguess = [0;0;0;0];
ew = .00001;
ev = .00001;

Toba = [1 0 0 0;0 1 0 0;0 0 1 .5;0 0 0 1];

Tobb = [0 -1 0 0;1 0 0 1;0 0 1 0;0 0 0 1];

Tobc = [0 1 0 2;0 0 1 0;1 0 0 1;0 0 0 1];

Tobd = [sqrt(3)/2 -1/2 0 3/sqrt(2);1/2 sqrt(3)/2 0 3/sqrt(2);0 0 1 0;0 0 0 1];

[thetalist, success] = IKinSpace(S, M, Tobd, thetaguess, ew, ev)


