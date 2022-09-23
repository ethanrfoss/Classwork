
clear all; clc;
%% number 5
x = @(t) cos(2*pi*20*t)+1/9*cos(3*2*pi*20*t)+1/25*cos(5*2*pi*20*t);
tsample = [0:1/44100:1/5];
xdiscrete = x(tsample);
xdiscreterand = xdiscrete + randn(1,length(xdiscrete));
Xdiscreterand = fft(xdiscreterand);

figure(1); hold on;
subplot(2,1,1);
plot(tsample,xdiscrete); % normal signal
xlabel('Time(sec)');
ylabel('xn');
title('Original Signal');
subplot(2,1,2);
plot(tsample,xdiscreterand); % signal with noise
xlabel('Time(sec)');
ylabel('xn');
title('Signal with Noise');

figure(2); hold on;
title('Discrete Fourier Transform of Noisy signal and Threshold');
plot(0:8820,abs(Xdiscreterand)/8820);
xlabel('Frequency(Hz)');
ylabel('Xn magnitude');
threshold = .035; %threshold determined by analyzing graph
plot(0:8820,threshold*ones(1,length(tsample)));

for i = 1:length(Xdiscreterand)
    if Xdiscreterand(i)/8820<threshold
        Xdiscreterand(i) = 0;
    end
end
xnoiseremoved = ifft(Xdiscreterand);

figure(3); hold on;
subplot(2,1,1);
plot(tsample,xdiscreterand);
title('Noisy Signal');
axis([0 0.2 -5 5]);
subplot(2,1,2);
plot(tsample, xnoiseremoved);
title('Signal with noise manually removed');



