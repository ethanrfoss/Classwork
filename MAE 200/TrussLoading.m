function [xe,fe] = TrussLoading(g,z)


[K,C,D] = TrussK(g);

[A,b,p] = GaussianElimination(K,z)

xe(10) = b(10)/A(10,10);
for i = 9:-1:1
    xe(i) = (b(i) - sum((xe(i+1:end).*A(i,i+1:end))))/A(i,i);
end

xe = xe';

fe = -C*D*xe;