function [U,SO,V,round] = sysjac(B,U,sort)
% SVD decomposition based on systolic jacobi iteration algorithm
% using round-robin array
% round:return the all jacobi iteration round 
[a,b]=size(B);
S=B;
%U=eye(a);
SO=eye(a);
V=eye(b);
TOL=1.e-4;
converge =TOL+1;
Q=[2:2:a];
P=[1:2:a-1];
round=0;
while converge>TOL
    converge =0;
    for i=1:a-1
    for n=1:a/2
    [S,U,error]=jacobiunit(P(n),Q(n),S,U,sort);
    converge=max(converge,error);
    end
    %decide next P,Q- round-robin array
    IQ=cat(2,Q(2:a/2),P(a/2));
    P=cat(2,1,Q(1),P(2:a/2-1));
    Q=IQ;
    end
    round=round+1;
    if(converge<TOL)
        break
    end
end
U=U';
%dervie S and V only for verification
for i = 1:a
    SO(i,i) = norm(S(i,:));
    V(i,:) = S(i,:)/SO(i,i);
end

end

