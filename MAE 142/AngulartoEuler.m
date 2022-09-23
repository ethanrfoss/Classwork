function Tw2e = AngulartoEuler(psi,theta,phi)

Tw2e = [0,sin(phi)*sec(theta),cos(phi)*sec(phi);0,cos(phi),-sin(phi);1 sin(phi)*tan(theta),cos(phi)*tan(theta)];

end
