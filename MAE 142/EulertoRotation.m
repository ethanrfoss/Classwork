
function TL2B = EulertoRotation(psi,theta,phi)

TL2B =  [cos(psi)*cos(theta),cos(theta)*sin(psi),-sin(theta); cos(psi)*sin(theta)*sin(phi)-cos(phi)*sin(psi),cos(psi)*cos(phi)+sin(psi)*sin(theta)*sin(phi),cos(theta)*sin(phi); sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta),cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi),cos(theta)*cos(phi)];

end