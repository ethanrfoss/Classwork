%% HW 5

%% Problem 5
load('MAE108S22p5data.mat');

mu = mean(eta);
v = var(eta);

figure; hold on; title('Problem 5: mean = -2.878E-19 | std = .4264');
histogram(eta,'Normalization','pdf');
plot(linspace(min(eta),max(eta),1000),normpdf(linspace(min(eta),max(eta),1000),mu,sqrt(v)));

%% Problem 6

alpha = 0;
for n = 1:length(zerouph)
    alpha = alpha + zerouph(n)^2;
end
alpha = sqrt(alpha/n)

figure; hold on; title('Problem 6: \alpha = 1.1485');
histogram(zerouph,'Normalization','pdf')
h = linspace(min(zerouph),max(zerouph),1000);
plot(h,2*h/alpha^2.*exp(-(h/alpha).^2));