figure(1); title('a');

bode(tf([1],[1 1]));

figure(2); title('b');

bode(tf([1 -1],[1]));

figure(3); title('c');

bode(tf([1 -1],[1 1]));

figure(4); title('d');

bode(tf([1 101 100],[1 0]));

figure(5); title('e');

bode(tf([1 0],[1 101 100]));

figure(6); title('f');

bode(tf([1],[1 1 1]));