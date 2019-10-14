function [U,SO,V,round] = ringjacob(B,U,sort)
% SVD decomposition based on ring jacobi iteration algorithm
% using ring-array
% round:return the all jacobi iteration round 
[a,b]=size(B);
S=B;
SO=eye(a);
V=eye(b);
TOL=1.e-4;
converge =TOL+1;
Q=[2:2:a];
P=[1:2:a-1];
round=0;
pol=1;  %polarity of ring order
while converge>TOL
    converge =0;
    for i=1:a-1
    for n=1:a/2
    [S,U,error]=jacobiunit(P(n),Q(n),S,U,sort);
    converge=max(converge,error);
    end
    %decide next ring array
    [P,Q]=nextpq(P,Q,i,pol);
    end
    round=round+1;
    if(pol==1)%change the polarity
        pol=0;
    else
        pol=1;
    end
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

