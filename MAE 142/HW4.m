%% HW 4

%% 4.2
G = tf([.25*42],[1 2.4 0]);

bode(G);

40.3*pi/180/2.83

%% 4.3
a = 20;
N = @(A) 2/pi*(asin(a/A)+(a/A)*sqrt(1-(a/A)^2));
K = .2;
G = 11*tf([1 .76],[1 1.47 -.32 0]);
T = @(A) minreal(G*N(A)/(1+G*K*N(A)));

Av = [80:.01:100];
for i = 1:length(Av)
    Av(i)
    if sum(real(pole(T(Av(i))))>=0)
        Av(i)
        pole(T(Av(i)))
        break;
    end
end

%% 4.5
syms Kq
solve(2*.08==(1.73+2.03*Kq)/sqrt(2.38+1.3195*Kq))

G = -2.03*tf([1 .65],[1 1.73 2.38]);
T = @(K) minreal(G/(1-G*K));
damp(T(.68))
