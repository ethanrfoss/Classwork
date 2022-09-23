

Rsc = [0 1 0; -1 0 0; 0 0 1];
Rcb = [0 -sqrt(3)/2 1/2; 1 0 0; 0 1/2 sqrt(3)/2];
Rsd = [2/3 -1/3 2/3;2/3 2/3 -1/3;-1/3 2/3 2/3];

Rdb = Rsd^-1*Rsc*Rcb;

theta = acos(.5*(trace(Rdb)-1));

wh = 1/(2*sin(theta))*(Rdb-Rdb^-1);
w =[-wh(2,3);wh(1,3);-wh(1,2)]


Rsb = Rsc*Rcb;
figure; hold on;
quiver3(0,0,0,1,0,0,1,'k')
quiver3(0,0,0,0,1,0,1,'k')
quiver3(0,0,0,0,0,1,1,'k')

quiver3(0,0,0,Rsb(1,1),Rsb(2,1),Rsb(3,1),1,'b');
quiver3(0,0,0,Rsb(1,2),Rsb(2,2),Rsb(3,2),1,'b');
%quiver3(0,0,0,Rsb(1,3),Rsb(2,3),Rsb(3,3),1,'b');

quiver3(0,0,0,Rsd(1,1),Rsd(2,1),Rsd(3,1),1,'r');
quiver3(0,0,0,Rsd(1,2),Rsd(2,2),Rsd(3,2),1,'r');
%quiver3(0,0,0,Rsd(1,3),Rsd(2,3),Rsd(3,3),1,'r');

xlabel('X');
ylabel('Y');
zlabel('Z');



figure; hold on;
quiver3(0,0,0,1,0,0,1,'b')
quiver3(0,0,0,0,1,0,1,'b')
quiver3(0,0,0,0,0,1,1,'b')

quiver3(0,0,0,Rdb(1,1),Rdb(2,1),Rdb(3,1),1,'k');
quiver3(0,0,0,Rdb(1,2),Rdb(2,2),Rdb(3,2),1,'r');
quiver3(0,0,0,Rdb(1,3),Rdb(2,3),Rdb(3,3),1,'r');

quiver3(0,0,0,w(1),w(2),w(3),'y');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;

dot(Rsb(1,:),w)
dot(Rsd(1,:),w)





v = [1;0;0];
vrot = v*cos(theta)+cross(w,v)*sin(theta)+w*(dot(w,v))*(1-cos(theta));