function [A,b,p] = GaussianElimination(A,b)

[m,n] = size(A);

if m~=n || m~=length(b)
    disp('A is not square, exiting');
    return;
end

%% Organize so that zeros are on bottom:
p = 1:m';
for i = 1:m-1
    [M,I] = sort(abs(A(i:m,i)),'descend');
    A(i:m,:) = A(I+i-1,:);
    b(i:m) = b(I+i-1);
    p(i:m) = p(I+i-1);
end

%% Row Reduce:
for i = 1:m-1
    for j = i+1:m
        t = A(j,i)/A(i,i);
        A(j,:) = A(j,:) - t*A(i,:);
        b(j) = b(j) - t*b(i);
    end
end

%% Make Ones:

for i = 1:m
    t = A(i,i);
    A(i,:) = A(i,:)/t;
    b(i) = b(i)/t;
end

end