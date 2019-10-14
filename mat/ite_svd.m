function [IU] = ite_svd(A,leftnumber,ite_time)
%iterative way to solve left singular vector-only matrix multiplication
P0=A*A';
PS=size(P0,1);
%%iterative way to solve singular vector one by one
IU=zeros(PS,leftnumber);

for a=1:leftnumber
    P=P0;
for i=1:ite_time
    P=(P*P)/norm(max(P));
end

IU(:,a)=P(:,1)/norm(P(:,1));
P0=(eye(PS)-IU(:,a)*IU(:,a)')*P0;
end

end


%{
error=(norm(IU)-norm(U))/norm(U)
errorm=abs(IU)-abs(U);
errorm=abs(errorm)./abs(U);
%}
%round
%{
A=rand(6,6);
[U,S,V]=svd(A);
UI=cell(6,1);%U's column vector
SI=zeros(6,1);
R=zeros(6,6);
for i=1:6
    UI{i}=U(:,i);
    SI(i)=S(i,i);
    R=R+SI(i)^16*UI{i}*UI{i}';
end
P=A*A';

for i=1:3
    P=P*P;
end
%test
RT=zeros(6,6);
for i=1:6
    RT(:,i)=SI(i)*U(:,i);
end
UT=zeros(6,6);
UT(1,:)=UI{1}';
ALL=RT*U';
SINGLE=RT*UT;
%}




