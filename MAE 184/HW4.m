%% HW 4 MAE 184

%% P1:

V = [133 137 141 143 145 147 149 151];
p230 = [109.68 119.44 129.779 135.168 140.707 146.398 152.243 158.244]*100/230;
fc = [10.37 11.15 11.97 12.40 12.84 13.29 13.74 14.21];

manp230 = [53 57 61 63 66 68 70 72];
mafc = [9.8 10.5 11 11.5 11.7 12.2 12.6 13];

figure; hold on;
plot(V,p230);
plot(V,manp230);
legend('JSBSim','Manufacturer Data');
title('Brake Horsepower vs. Airspeed');
xlabel('True Airspeed[kts]'); ylabel('Brake Horsepower[Percentage of 230 hp]');

figure; hold on;
plot(V,fc);
plot(V,mafc);
legend('JSBSim','Manufacturer Data');
title('Fuel Consumption Rate vs. Airspeed');
xlabel('True Airspeed[kts]'); ylabel('Fuel Consumption Rate[gal/hr]');

%% P2:
g = 32.2;
TSL = 518.7;
a = -.00356;
R = 1716;

h = [500 1500 2500 3500 4500 5500];
vc = 9620.52;

r = (1+a*h/TSL).^(-(1+g/(a*R)));

vt = vc./sqrt(r);
gamma = asin(680./vt);

fc = [10.83222 10.904981 10.980915 11.059855 11.141892 11.227434];

t = 1000/680;
f = zeros(6,1);
for i = 1:length(fc)
    f(i+1) = f(i)+fc(i)*t/60;
end

manf = [0 .5 1 1.5 2 2.6 3.1];
h = (0:6)*1000;

figure; hold on;
plot(f,h);
plot(manf,h);
legend('JSBSim','Manufacturer Data');
title('Altitude vs. Fuel Consumption');
xlabel('Fuel Used[Gallons]');
ylabel('Altitude[ft]');




