A = csvread('compass.csv');

t = A(:,1);
ax = A(:,2);
ay = A(:,3);
az = A(:,4);
mx = A(:,5);
my = A(:,6);
mz = A(:,7);

for i = 1:length(t)
    [psi(i),theta(i),phi(i)] = AccelandMagtoEuler(ax(i),ay(i),az(i),mx(i),my(i),mz(i));
end

figure; hold on;

plot(t,psi*180/pi,'r');
plot(t,theta*180/pi,'b');
plot(t,phi*180/pi,'g');

legend('Heading(\psi)','Pitch(\theta)','Roll(\phi)');
title('Euler Angles')
xlabel('Time[sec]');
ylabel('Euler Angles[deg]');