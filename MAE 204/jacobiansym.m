function Js = jacobiansym(S)

[m,n] = size(S);
if m~=6
    disp('Incorrect S dimensions');
    return;
end

Etot = sym(eye(4,4));
for i = 1:n
    
    eval(sprintf('syms theta%1d',i))
    theta(i) = eval(sprintf('theta%1d',i));
    
    Js(:,i) = Adjoint(Etot)*S(:,i);
    
    w = [0 -S(3,i) S(2,i); S(3,i) 0 -S(1,i); -S(2,i) S(1,i) 0];
    topleft = eye(3,3) + sin(theta(i))*w+(1-cos(theta(i)))*w*w;
    topright = (eye(3,3)*theta(i)+(1-cos(theta(i)))*w+(theta(i)-sin(theta(i)))*w*w)*S(4:6,i);
    matexp = [topleft,topright;zeros(1,3),1];

    Etot = Etot*matexp;
    

end

Js = simplify(Js);