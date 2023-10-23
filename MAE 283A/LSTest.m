num = 10^3*[.3855 -1.3785 1.9088 -1.1858 .2476 .0257];
den = [1 -5.3705 12.4447 -15.9157 11.8458 -4.8652 .8618];
N = 1000; Ts = .9;
G0 = tf(num,den,1);
w = linspace(0,pi,N)';
[mag,phase] = bode(G0,w);

Gspa = squeeze(mag).*exp(squeeze(phase)*pi/180*1i);

X=[exp(-[length(den)-length(num):length(den)-1]*1j.*w) -exp(-[1:length(den)-1]*1j.*w).*Gspa];
theta=[real(X);imag(X)]\[real(Gspa);imag(Gspa)];

Ghat = tf(theta(1:length(num))',[1 theta(length(num)+1:end)'],Ts);
[mhat,phat] = bode(Ghat,w);

figure(3);
subplot(2,1,1);
loglog(w,squeeze(mag),'g',w,squeeze(mhat),'r'); legend('G_{spa}','Ghat');
ylabel('Magnitude');
subplot(2,1,2);
semilogx(w,squeeze(phase),'g',w,squeeze(phat),'r');
xlabel('Frequency[rad/s]'); ylabel('Phase');