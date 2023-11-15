
%% Midterm


%% Problem 5
f = 2^(1/3);

c1 = 1/(2*(2-f));
c4 = c1;
c2 = (1-f)*c1;
c3 = c2;
d1 = 1/(2-f);
d3 = d1;
d2 = -f*d1;

e11 = d1+d2+d3
e12 = c1+c2+c3+c4
e21 = d1*c1+d2*c1+d2*c2+d3*c1+d3*c2+d3*c3
e22 = d1*c2+d1*c3+d2*c3+d1*c4+d2*c4+d3*c4
e31 = c2*d1*c1+c3*d1*c1+c3*d2*c1+c3*d2*c2+c4*d1*c1+c4*d2*c1+c4*d2*c2+c4*d3*c1+c4*d3*c2+c4*d3*c3
e32 = d2*c2*d1+d3*c2*d1+d3*c3*d1+d3*c3*d2
e41 = d2*c2*d1*c1+d3*c2*d1*c1+d3*c3*d1*c1+d3*c3*d2*c1+d3*c3*d2*c2
e42 = d2*c3*d1*c2+d2*c4*d1*c2+d3*c4*d1*c2+d3*c4*d1*c3+d3*c4*d2*c3

%% Problem 6

% SI1:

syms h k

c1 = 1;
d1 = 1;

D1a = (eye(2,2) + h*c1*[0 1;0 0]);
D1b = (eye(2,2) + h*d1*[0 0;-k 0]);

D1 = D1b*D1a;

sigma = [1/2*((D1(1,1)+D1(2,2)+sqrt(4*D1(2,1)*D1(1,2)+(D1(1,1)-D1(2,2))^2)));1/2*((D1(1,1)+D1(2,2)-sqrt(4*D1(2,1)*D1(1,2)+(D1(1,1)-D1(2,2))^2)))];
sigma = matlabFunction(sigma,'Vars',{h,k});

NP = 201;
B = [-3;3;-3;3];
LR = linspace(B(1),B(2),NP);
LI = linspace(B(3),B(4),NP);

for i = 1:NP
    for j = 1:NP
        L = LR(i)+1i*LI(j);
        sig(j,i) = max(abs(sigma(L,1)));
    end
    i
end
figure; contourf(LR,LI,1./sig,[.9999999999 1.000000000001],'k-'), colormap autumn, axis('square'), hold on
%plot([B(1) B(2)],[0,0],'k-'), plot([0,0],[B(3) B(4)],'k-');
title('Stability contours for SI1');
xlabel('L (real)')
ylabel('L (imaginary)')

% SI2:

syms h k

c1 = 1/2;
c2 = 1/2;
d1 = 1;

D1a = (eye(2,2) + h*c1*[0 1;0 0]);
D1b = (eye(2,2) + h*d1*[0 0;-k 0]);
D2a = (eye(2,2) + h*c2*[0 1;0 0]);

D2 = D2a*D1b*D1a;

sigma = [1/2*((D2(1,1)+D2(2,2)+sqrt(4*D2(2,1)*D2(1,2)+(D2(1,1)-D2(2,2))^2)));1/2*((D2(1,1)+D2(2,2)-sqrt(4*D2(2,1)*D2(1,2)+(D2(1,1)-D2(2,2))^2)))];
sigma = matlabFunction(sigma,'Vars',{h,k});

NP = 201;
B = [-3;3;-3;3];
LR = linspace(B(1),B(2),NP);
LI = linspace(B(3),B(4),NP);

for i = 1:NP
    for j = 1:NP
        L = LR(i)+1i*LI(j);
        sig(j,i) = max(abs(sigma(L,1)));
    end
    i
end
figure; contourf(LR,LI,1./sig,[.9999999999 1.000000000001],'k-'), colormap autumn, axis('square'), hold on
%plot([B(1) B(2)],[0,0],'k-'), plot([0,0],[B(3) B(4)],'k-');
title('Stability contours for SI2');
xlabel('L (real)')
ylabel('L (imaginary)')

% SI4:

syms h k

f = 2^(1/3);

c1 = 1/(2*(2-f));
c4 = c1;
c2 = (1-f)*c1;
c3 = c2;
d1 = 1/(2-f);
d3 = d1;
d2 = -f*d1;

D1a = (eye(2,2) + h*c1*[0 1;0 0]);
D1b = (eye(2,2) + h*d1*[0 0;-k 0]);
D2a = (eye(2,2) + h*c2*[0 1;0 0]);
D2b = (eye(2,2) + h*d2*[0 0;-k 0]);
D3a = (eye(2,2) + h*c3*[0 1;0 0]);
D3b = (eye(2,2) + h*d3*[0 0;-k 0]);
D4a = (eye(2,2) + h*c4*[0 1;0 0]);

D4 = D4a*D3b*D3a*D2b*D2a*D1b*D1a;

sigma = [1/2*((D4(1,1)+D4(2,2)+sqrt(4*D4(2,1)*D4(1,2)+(D4(1,1)-D4(2,2))^2)));1/2*((D4(1,1)+D4(2,2)-sqrt(4*D4(2,1)*D4(1,2)+(D4(1,1)-D4(2,2))^2)))];
sigma = matlabFunction(sigma,'Vars',{h,k});

NP = 201;
B = [-3;3;-3;3];
LR = linspace(B(1),B(2),NP);
LI = linspace(B(3),B(4),NP);

for i = 1:NP
    for j = 1:NP
        L = LR(i)+1i*LI(j);
        sig(j,i) = max(abs(sigma(L,1)));
    end
    i
end
figure; contourf(LR,LI,1./sig,[.9999999999 1.000000000001],'k-'), colormap autumn, axis('square'), hold on
%plot([B(1) B(2)],[0,0],'k-'), plot([0,0],[B(3) B(4)],'k-');
title('Stability contours for SI4');
xlabel('L (real)')
ylabel('L (imaginary)')


