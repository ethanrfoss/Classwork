
P = tf([1],[1]);
I = tf([300],[1 0]);
D = tf([50 0],[1]);
PC = tf([1],[1 -20]);
%% 6.8a
rlocus(tf([1],[1 3 2]),tf([0]))


%% 6.8b
G = tf([1 -1],[1 1]);
rlocus(G*(D));
sisotool(G)

%% 6.8C

G = tf([1],[1 -3 2]);
rlocus(G*(P+D))

%% 6.8d

G = tf([1 -1],[1 -2])
rlocus(G*(PC),20:30)

%% 6.8g

G = tf([1],[1 2 1 2])
rlocus(G*(D))

%% 6.15

G = tf([1/1840],[1 240/1840]);
K = tf([1 1],[1 0]);
rlocus(G*K);