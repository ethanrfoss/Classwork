function R=RHS(x,u,s)
E=ComputeE(x,s);
N=ComputeN(x,u,s);
R=E\N;
end