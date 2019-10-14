%this file is for the ring/systolic jacobi array comparision, 
RS=zeros(1,256-32);
RD=zeros(1,256-32);
RATE=zeros(1,128);
cdim=32;

for odim=32:16:320 % test different matrix size,suppose the core dimension are 32
A=rand(odim,cdim*cdim);
UI=eye(odim);
[U,S,V]=svd(A);
[a,~]=size(S);

%% systolic jacobi iteration way
[UO,SO,VO,rounds] = sysjac(A,UI,1);
RS(odim)=rounds;
errors=(norm(S(1:a,1:a)-SO(1:a,1:a))^2/norm(S)^2)*100;

%% ring jacobi iteration way
[UR,SR,VR,roundr] = ringjacob(A,UI,1);
RD(odim)=roundr;
errorr=(norm(SR(1:a,1:a)-S(1:a,1:a))^2/norm(S)^2)*100;

%% direct iteration way
[IR]=ite_svd(A,cdim,15);
errori=norm(abs(IR)-abs(U(:,1:cdim)))^2/norm(U)^2;
errorc=(abs(IR)-abs(U(:,1:cdim)))/abs(U(:,1:cdim));

%% print result
fprintf('matrix_size: %d sys_round = %d, error = %d %',odim,rounds,errors);
fprintf(',ring_round= %d ,error= %d\n',roundr,errorr);

end


