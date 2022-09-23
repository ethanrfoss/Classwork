function [T] = LatLongTransformation(lat,long)

T = [-sind(lat)*cosd(long),-sind(long),-cosd(lat)*cosd(long);-sind(lat)*sind(long),cosd(long),-cosd(lat)*sind(long);cosd(lat),0,-sind(lat)];

end