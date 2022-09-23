function R=RHSX(x,s,X,Q,R)
E=ComputeE(x,s);
A=ComputeA(x,s);
R = -(E')^-1*A'*X-X*A*E^-1+X*s.B*R^-1*s.B'*X-(E')^-1*Q*E^-1;
end