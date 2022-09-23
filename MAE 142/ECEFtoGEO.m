function [lat, long, h] = ECEFtoGEO(X,Y,Z)

lat = atan2(Z,sqrt(X^2+Y^2))*180/pi;
long = atan2(Y,X)*180/pi;
h = sqrt(X^2+Y^2+Z^2) - 6367*10^3;


end