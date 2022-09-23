function [PECEF,PGEO,b] = GPS(s1,p1,s2,p2,s3,p3,s4,p4)

alpha = [s1'*s1-p1^2; s2'*s2-p2^2; s3'*s3-p3^2; s4'*s4-p4^2];
H = [2*s1' , 2*p1; 2*s2' , 2*p2; 2*s3' , 2*p3; 2*s4' , 2*p4];
e = [1;1;1;1];

Hinv = H^-1;
P = Hinv(1:3,:);
Q = Hinv(4,:);

a = e'*(P'*P-Q'*Q)*e;
b = 2*alpha'*(P'*P-Q'*Q)*e-1;
c = alpha'*(P'*P-Q'*Q)*alpha;

beta = roots([a b c]);

XB1 = Hinv*(alpha+e*beta(1));
XB2 = Hinv*(alpha+e*beta(2));

P = [XB1(1),XB2(1);XB1(2),XB2(2);XB1(3),XB2(3)];
b = [XB1(4),XB2(4)];

[lat1,long1,h1] = ECEFtoGEO(P(1,1),P(2,1),P(3,1)); 
[lat2,long2,h2] = ECEFtoGEO(P(1,2),P(2,2),P(3,2)); 

if abs(h1)<abs(h2)
    PECEF = P(:,1);
    PGEO = [lat1,long1,h1]';
    b = b(1);
else
    PECEF = P(:,2);
    PGEO = [lat2,long2,h2]';
    b = b(2);
end