
tfb = (0:.1:20)/1000;
r3 = .5/100 + tfb;
r2 = .5/100;
r1 = .4/100;
hi = 100;
ho = 6;
ka = 25;
kfg = .06;
for i = 1:length(tfb)
q(i) = 75/(1/(2*pi*r1*hi)+log(r2/r1)/(2*pi*ka)+log(r3(i)/r2)/(2*pi*kfg)+1/(2*pi*r3(i)*ho));
end

figure(1);
plot(tfb,q);

xlabel('Fiber Glass Thickness(m)');
ylabel('Heat Loss per Unit Length(W/m)');