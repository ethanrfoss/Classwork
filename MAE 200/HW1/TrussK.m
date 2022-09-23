function [K,C,D] = TrussK(g)

xa = 0; ya = 0;
xb = 0; yb = 1;
xg = 3; yg = 1;
c = 10^7;

l = [sqrt((g(2)-ya)^2+(g(1)-xa)^2);sqrt((g(4)-ya)^2+(g(3)-xa)^2);sqrt((g(4)-yb)^2+(g(3)-xb)^2);sqrt((g(4)-g(2))^2+(g(3)-g(1))^2);...
    sqrt((g(6)-g(2))^2+(g(5)-g(1))^2); sqrt((g(8)-g(2))^2+(g(7)-g(1))^2); sqrt((g(8)-g(4))^2+(g(7)-g(3))^2); sqrt((g(8)-g(6))^2+(g(7)-g(5))^2);...
    sqrt((yg-g(6))^2+(xg-g(5))^2); sqrt((yg-g(8))^2+(xg-g(7))^2)];

C = diag([c/l(1) c/l(2) c/l(3) c/l(4) c/l(5) c/l(6) c/l(7) c/l(8) c/l(9) c/l(10)]);

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


D = [cos(t1) sin(t1) 0 0 0 0 0 0 0 0;...
    0 0 cos(t2) sin(t2) 0 0 0 0 0 0;....
    0 0 cos(t3) sin(t3) 0 0 0 0 0 0;...
    -cos(t4) -sin(t4) cos(t4) sin(t4) 0 0 0 0 0 0;...
    -cos(t5) -sin(t5) 0 0 cos(t5) sin(t5) 0 0 0 0;...
    -cos(t6) -sin(t6) 0 0 0 0 cos(t6) sin(t6) 0 0;...
    0 0 -cos(t7) -sin(t7) 0 0 cos(t7) sin(t7) 0 0;...
    0 0 0 0 -cos(t8) -sin(t8) cos(t8) sin(t8) 0 0;...
    0 0 0 0 -cos(t9) -sin(t9) 0 0 cos(t9) sin(t9);...
    0 0 0 0 0 0 -cos(t10) -sin(t10) cos(t10) sin(t10)];

K = D'*C*D;

end
