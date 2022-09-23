G1 = tf([2],[1 2]);
G2 = tf([0.5],[1 0]);
G3 = tf([5],[1 0]);
G = G1*G2*G3;
K = tf([.3 .13],[1]);

S = minreal(1/(1+G*K));
S0 = evalfr(S,0)

D = minreal(G/(1+G*K))
D0 = evalfr(D,0)

error = 1*S0 + .5*D0;

%% Sin

[magy, phasey] = bode(S,1)
[magv, phasev] = bode(D,2)
