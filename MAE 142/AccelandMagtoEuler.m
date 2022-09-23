function [psi,theta,phi] = AccelandMagtoEuler(ax,ay,az,mx,my,mz)

theta = atan2(ax,sqrt(ay^2+az^2));
phi = atan2(-ay,-az);

dec = 12*pi/180;
psi = dec+atan2(mz*sin(phi)-my*cos(phi),mx*cos(theta)+my*sin(theta)*sin(phi)+mz*sin(theta)*cos(phi));
end