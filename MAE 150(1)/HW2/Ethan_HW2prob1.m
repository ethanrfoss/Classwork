%%
clear all; clc;
mew = zeros(1,5000);
mewsin = zeros(1,5000);
n=1:5000;
for i = 1:length(n)
    for j=1:i
        mew(i) = mew(i) + (rand)^3/n(i);
        mewsin(i) = mewsin(i) + sin(pi*rand)/n(i);
    end
end
figure(1); hold on;
subplot(2,1,1); hold on;
plot(n,mew);
plot(n, .25*ones(1,length(n)));
legend('Estimator','Theoretical');
title('Monte Carlo Estimator for x^3');
xlabel('n'); ylabel('mew');

subplot(2,1,2); hold on;
plot(n,mewsin);
plot(n, 2/pi*ones(1,length(n)));
axis([0 5000 0 1]);
legend('Estimator','Theoretical');
title('Monte Carlo Estimator for sin(pi*x)');
xlabel('n'); ylabel('mew');

%% If the random integral values are interdependent
%clear all; clc;
%mewtot = 0;
% mewsintot = 0;
% mew = zeros(1,5000);
% mewsin = zeros(1,5000);
% n=1:5000;
% for i = 1:length(n)
%     mewtot = mewtot + (rand)^3;
%     mew(i) = mewtot/n(i);
%     mewsintot = mewsintot + sin(pi*rand);
%     mewsin(i) = mewsintot/n(i);
% end
% figure(1); hold on;
% subplot(2,1,1); hold on;
% plot(n,mew);
% plot(n, .25*ones(1,length(n)));
% legend('Estimator','Theoretical');
% title('Monte Carlo Estimator for x^3');
% xlabel('n'); ylabel('mew');
% 
% subplot(2,1,2); hold on;
% plot(n,mewsin);
% plot(n, 2/pi*ones(1,length(n)));
% axis([0 5000 0 1]);
% legend('Estimator','Theoretical');
% title('Monte Carlo Estimator for sin(pi*x)');
% xlabel('n'); ylabel('mew');