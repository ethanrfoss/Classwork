function m = MagnetictoLocal(dip,dec,M)

m = [cos(dec),-sin(dec),0;sin(dec),cos(dec),0;0,0,1]*[cos(dip),0,-sin(dip);0,1,0;sin(dip),0,cos(dip)]*[M;0;0];

end