function T = ArmTransformBody(M,B,theta)

[~,n] = size(B);

if n~=length(theta)
    disp('Number of Inputs not equal to number of joints');
    return;
end

T = M;

for i = 1:n
    w = [0 -B(3,i) B(2,i); B(3,i) 0 -B(1,i); -B(2,i) B(1,i) 0];
    topleft = eye(3,3) + sin(theta(i))*w+(1-cos(theta(i)))*w*w;
    topright = (eye(3,3)*theta(i)+(1-cos(theta(i)))*w+(theta(i)-sin(theta(i)))*w*w)*B(4:6,i);
    matexp = [topleft,topright;zeros(1,3),1];
    
    T = T*matexp;
    
end

end