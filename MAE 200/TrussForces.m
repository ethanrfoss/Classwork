function f = TrussForces(g)

xa = 0; ya = 0;
xb = 0; yb = 1;
xg = 3; yg = 1;
W = 1000;

t1 = atan2(g(2)-ya,g(1)-xa);
t2 = atan2(g(4)-ya,g(3)-xa);
t3 = atan2(g(4)-yb,g(3)-xb);
t4 = atan2(g(4)-g(2),g(3)-g(1));
t5 = atan2(g(6)-g(2),g(5)-g(1));
t6 = atan2(g(8)-g(2),g(7)-g(1));
t7 = atan2(g(8)-g(4),g(7)-g(3));
t8 = atan2(g(8)-g(6),g(7)-g(5));
t9 = atan2(yg-g(6),xg-g(5));
t10 = atan2(yg-g(8),xg-g(7));

A = [cos(t1) 0 0 -cos(t4) -cos(t5) -cos(t6) 0 0 0 0;...
     sin(t1) 0 0 -sin(t4) -sin(t5) -sin(t6) 0 0 0 0;...
     0 cos(t2) cos(t3) cos(t4) 0 0 -cos(t7) 0 0 0;...
     0 sin(t2) sin(t3) sin(t4) 0 0 -sin(t7) 0 0 0;...
     0 0 0 0 cos(t5) 0 0 -cos(t8) -cos(t9) 0;...
     0 0 0 0 sin(t5) 0 0 -sin(t8) -sin(t9) 0;...
     0 0 0 0 0 cos(t6) cos(t7) cos(t8) 0 -cos(t10);...
     0 0 0 0 0 sin(t6) sin(t7) sin(t8) 0 -sin(t10);...
     0 0 0 0 0 0 0 0 cos(t9) cos(t10);...
     0 0 0 0 0 0 0 0 sin(t9) sin(t10)];
 
 b = [0;0;0;0;0;0;0;0;0;W];
 
 [A,b,p] = GaussianElimination(A,b);
 
f(10) = b(10)/A(10,10);
for i = 9:-1:1
    f(i) = (b(i) - sum((f(i+1:end).*A(i,i+1:end))))/A(i,i);
end


end