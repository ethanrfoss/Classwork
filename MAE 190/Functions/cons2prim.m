function [rho,u,v,T,p,e,Et] = cons2prim(U,R,cv)

rho = squeeze(U(1,:,:));
Et = squeeze(U(4,:,:));
u = squeeze(U(2,:,:))./rho;
v = squeeze(U(3,:,:))./rho;
e = Et./rho-(u.^2+v.^2)/2;
T = e./cv;
p = rho.*R.*T;

end