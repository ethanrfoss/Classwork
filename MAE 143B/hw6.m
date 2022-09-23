syms s
Y = -5/100*exp(-.000001*s)/s+5/100*exp(-.036*s)/s

G = -s^2/(s^2+2.51327*s+246.74);

Z = G*Y;

fplot(ilaplace(Z));
title('Speed: 100 km/hr');
xlabel('time(s)');
ylabel('Z position(m)')
axis([0 2.5 -.05 .05]);