%% HW 6 Problem 3

S = [40 9 100 50 15 65 25 60 95 65 30 125];
SD = [15 2 40 15 4 25 5 25 30 24 8 45];

ba = regress(SD',[ones(length(S'),1) S']);

b = regress(SD',S');

figure; hold on;
plot(S,SD,'o');
plot(S,ba(1)+ba(2)*S);
plot(S,b*S);
legend('Actual','Regression With Intercept','Regression Without an Intercept');
xlabel('Speed'); ylabel('Stopping Distance'); title('Problem 3');

beta = sum((S-mean(S)).*(SD-mean(SD)))/sum((S-mean(S)).^2);
alpha = mean(SD)-beta*mean(S);

corcoef = sum((S-mean(S)).*(SD-mean(SD)))/sqrt(sum((S-mean(S)).^2)*sum((SD-mean(SD)).^2));

%% HW6 Problem 4

GNP = [600 2700 2900 4200 3100 5400 8600 10300];
E = [1000 700 1400 2000 2500 2700 2500 4000];
n = length(E);
[GNP,sl] = sort(GNP);
E = E(sl);

ba = regress(E',[ones(length(GNP'),1) GNP']);

syx = sqrt(1/(n-2)*(sum((E-mean(E)).^2)-ba(2)^2*sum((GNP-mean(GNP)).^2)));
up = ba(1)+ba(2)*GNP + tinv(.975,n-2)*syx*sqrt(1/n+(GNP-mean(GNP)).^2/sum((GNP-mean(GNP)).^2));
down = ba(1)+ba(2)*GNP - tinv(.975,n-2)*syx*sqrt(1/n+(GNP-mean(GNP)).^2/sum((GNP-mean(GNP)).^2));

figure; hold on;
plot(GNP,E,'o');
plot([0,GNP],ba(1) +ba(2)*[0,GNP]);
plot(GNP,up);
plot(GNP,down);
legend('Raw Data','Best Fit Line','Upper 95% Confidence Level','Lower 95% Confidence Level');
xlabel('GNP'); ylabel('Energy Consumption'); title('Problem 4');

corcoef = sum((E-mean(E)).*(GNP-mean(GNP)))/sqrt(sum((E-mean(E)).^2)*sum((GNP-mean(GNP)).^2));