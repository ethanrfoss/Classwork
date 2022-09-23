syms a e i n w E

da = 2*sin(i)/n^2*(sqrt(1-e^2)*cos(w)*sin(E)-sin(w)*(1-cos(E)));

int(da*(1-e*cos(E)),E,[0 2*pi])

de = sqrt(1-e^2)*sin(i)/(n^2*a)*(3/2*cos(w)*E-2*e*cos(w)*sin(E)+1/4*cos(w)*sin(2*E)-sqrt(1-e^2)/4*sin(w)*(1-cos(2*E)));

edot = 3/2*sqrt(1-e^2)/(n*a)*cos(w)*sin(i);

1/(2*pi)*int(de-1/n*edot*(E-e*sin(E))*(1-e*cos(E)),E,[0 2*pi])