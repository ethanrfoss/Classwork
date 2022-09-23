%%3 

long = atan2d(-19757,-1931)

%% 4 

b = (23582959-25081922)/299792458

%% 5

latSD = 32.734; longSD = -117.190;

latHON = 21.318; longHON = -157.920;

a = 6378;

[ESD] = GEOtoECEF(latSD,longSD,0)
[EHON] = GEOtoECEF(latHON,longHON,0)
ESD = ESD/(sqrt(sum(ESD.^2)))
EHON = EHON/(sqrt(sum(EHON.^2)))
dgc = a*acos(dot(ESD,EHON));

dgc = a*acos(sind(latSD)*sind(latHON) + cosd(latSD)*cosd(latHON)*cosd(longHON-longSD));


%% 6

t = 2460; 
d = 230*t;
latA = -1.36; longA = -32.83;
X = 27.8;

latB = 565800/(1000*secd(X)*(6378+11))+latA

TA = log(secd(latA)+tand(latA));
TB = log(secd(latB)+tand(latB))

longB = tand(X)*(TB-TA)+longA;

ECEFB = GEOtoECEF(latB,longA,11000);
ECEFA = GEOtoECEF(latA,longA,11000);

ECEFB = ECEFB/sqrt(sum(ECEFB.^2));
ECEFA = ECEFA/sqrt(sum(ECEFA.^2));

d = (a+11)*acos(dot(ECEFB,ECEFA))

%% 6.2

ECEFINTOL = GEOtoECEF(-1.36,-32.83,11000);

ECEFCRASH = ECEFINTOL + LatLongTransformation(-1.36,-32.83)*[230*2460*cosd(27.8);230*2460*sind(27.8); 0];

[latC, longC, hC] = ECEFtoGEO(ECEFCRASH(1),ECEFCRASH(2),ECEFCRASH(3));


a*acos(dot(ECEFINTOL,ECEFCRASH)/(sqrt(sum(ECEFINTOL.^2))*sqrt(sum(ECEFCRASH.^2))))


%% 7 
psi = 0; theta = 0; phi = 11;
TL2B = EulertoRotation(psi,theta,phi*pi/180);
pqr = [0;.7;2.9];
 p =pqr(1); q = pqr(2); r= pqr(3);

psthph = TL2B^-1*pqr

psidot = (q*sind(phi)+r*cosd(phi))/cosd(theta)
