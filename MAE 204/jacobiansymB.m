function Jb = jacobiansymB(B)

[m,n] = size(B);
if m~=6
    disp('Incorrect S dimensions');
    return;
end

Etot = sym(eye(4,4));
for i = n:-1:1
    
    eval(sprintf('syms theta%1d',i))
    theta(i) = eval(sprintf('theta%1d',i));
    
    Jb(:,i) = Adjoint(Etot)*B(:,i);
    
    w = -[0 -B(3,i) B(2,i); B(3,i) 0 -B(1,i); -B(2,i) B(1,i) 0];
    topleft = eye(3,3) + sin(theta(i))*w+(1-cos(theta(i)))*w*w;
    topright = (eye(3,3)*theta(i)+(1-cos(theta(i)))*w+(theta(i)-sin(theta(i)))*w*w)*-B(4:6,i);
    matexp = [topleft,topright;zeros(1,3),1];

    Etot = Etot*matexp;
    

end
Jb
Jb = simplify(Jb);

end