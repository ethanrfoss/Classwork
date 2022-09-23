%% HW1 MAE 108

%% 1:
figure(1);
xrandu = rand(1000,1);

subplot(2,2,1); plot(xrandu);

subplot(2,2,2); histogram(xrandu);
xlabel('Value'); ylabel('Occurences');

subplot(2,2,3); histogram(xrandu,[0:.25:1]);

subplot(2,2,4); histogram(xrandu,'Normalization','pdf');

figure(2);
xrandu = rand(10000,1);

subplot(2,2,1); plot(xrandu);

subplot(2,2,2); histogram(xrandu);
xlabel('Value'); ylabel('Occurences');

subplot(2,2,3); histogram(xrandu,[0:.25:1]);

subplot(2,2,4); histogram(xrandu,'Normalization','pdf');

figure(3);
xrandu = randn(1000,1);

subplot(2,2,1); plot(xrandu);

subplot(2,2,2); histogram(xrandu);
xlabel('Value'); ylabel('Occurences');

subplot(2,2,3); histogram(xrandu,[0:.25:1]);

subplot(2,2,4); histogram(xrandu,'Normalization','pdf');

figure(4);
xrandu = randn(10000,1);

subplot(2,2,1); plot(xrandu);

subplot(2,2,2); histogram(xrandu);
xlabel('Value'); ylabel('Occurences');

subplot(2,2,3); histogram(xrandu,[0:.25:1]);

subplot(2,2,4); histogram(xrandu,'Normalization','pdf');