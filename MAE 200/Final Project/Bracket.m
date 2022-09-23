function [AA,AB,AC,JA,JB,JC]= Bracket(ComputeJ,AA,AB,JA,X,P,V)%NumericalRenaissanceCodebase1.0
JB=ComputeJ(X+AB*P,V);
if JB>JA
    %[AA,AB]=Swap(AA,AB);
    temp = AA;
    AA = AB;
    AB = temp;
    %[JA,JB]=Swap(JA,JB);
    temp = JA;
    JA = JB;
    JB = temp;
end
AC=AB+2*(AB-AA);
JC=ComputeJ(X+AC*P,V);
end%functionBracket