d = randn(1,10000);
N = 10000;
u(1) = 4/10*d(1);
u(2) = 4/10*d(2)+3/10*d(1);
u(3) = 4/10*d(3)+3/10*d(2)+2/10*d(1);
for i = 4:length(d)
    u(i) = 4/10*d(i)+3/10*d(i-1)+2/10*d(i-2)+1/10*d(i-3);
end

tn = 10;
Rud = xcorr(u,d,tn,'biased');

figure(1); hold on; title('Cross Correlation R_u^N(T)');
plot([-tn:tn],Rud);
xlabel('\tau'); ylabel('R_u^N(\tau)');

UN = fft(u);

PUN = 1/N*UN.*conj(UN);

w = [0:.001:pi];
H = 4/10 + 3/10*exp(-1*1i*w) + 2/10*exp(-2*1i*w) + 1/10*exp(-3*1i*w);

figure(2);
subplot(2,1,1);
semilogx(w,abs(H)); title('Magnitude of H(z)'); 
subplot(2,1,2);
semilogx(w,angle(H)); title('Phase of H(z)'); xlabel('Frequency[rad/s]');

figure(3);hold on; title('Periodogram and Magnitude');
semilogx([1:N/2]*4*pi/(2*N+1),abs(PUN(1:end/2)));
semilogx(w,abs(H).^2); legend('P_u^N(w)','|H(jw)|^2');
xlabel('Frequency[rad/s]'); ylabel('Magnitude');





