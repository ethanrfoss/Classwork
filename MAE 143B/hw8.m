
z = 10;
p = 1;
k=1;

m=1;
b=5;

K = k*tf([1 z],[1 p]);
G = tf([1],[m b 0]);
Kpi = tf([10 1],[1 0]);
rlocus(G*Kpi);

bode(Kpi)
