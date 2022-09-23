function T = ArmTransformSpace(M,S,theta)

[~,n] = size(S);

if n~=length(theta)
    disp('Number of Inputs not equal to number of joints');
    return;
end

T = M;

for i = n:-1:1
    w = [0 -S(3,i) S(2,i); S(3,i) 0 -S(1,i); -S(2,i) S(1,i) 0];
    topleft = eye(3,3) + sin(theta(i))*w+(1-cos(theta(i)))*w*w;
    topright = (eye(3,3)*theta(i)+(1-cos(theta(i)))*w+(theta(i)-sin(theta(i)))*w*w)*S(4:6,i);
    matexp = [topleft,topright;zeros(1,3),1];
    
    T = matexp*T;
    
end