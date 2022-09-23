function [E] = GEOtoECEF(lat,long,h)

R = 6367*10^3;

X = (R+h)*cosd(lat)*cosd(long);
Y = (R+h)*cosd(lat)*sind(long);
Z = (R+h)*sind(lat);

E = [X;Y;Z];
end