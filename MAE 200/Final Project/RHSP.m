function R=RHSP(x,s,P,Q,R)
E=ComputeE(x,s);
A=ComputeA(x,s);
R = -(E')^-1*A'*P-P*A*E^-1+P*s.C'*R^-1*s.C*P-(E')^-1*Q*E^-1;
end