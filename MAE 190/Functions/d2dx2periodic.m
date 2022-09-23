function d2dx2 = d2dx2periodic(f,dx)

[nx,ny] = size(f);
d2dx2 = zeros(size(f));
for i = 2:nx-1
    for j = 1:ny
        d2dx2(i,j) = (f(i+1,j)-2*f(i,j)+f(i-1,j))/dx^2;
    end
end

i=1;
for j = 1:ny
    d2dx2(i,j) = (f(i+1,j)-2*f(i,j)+f(nx,j))/dx^2;
end

i=nx;
for j = 1:ny
    d2dx2(i,j) = (f(1,j)-2*f(i,j)+f(i-1,j))/dx^2;
end

end