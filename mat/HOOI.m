function [T]=HOOI(X,order,init,sort)
%HOOI with warm-start Jacobi itertaion,aims at 3-dim vector
normX=norm(X);
S=size(X);
xdim=ndims(X);
U=cell(xdim:1);

%initialization part
for i=1:3
    U{i}=eye(S(i));
    A(:,:,i)=init{i}(:,1:order);
end
tol=1e-4;
error=tol+1;
comptol=1e-4; %comparative error tolerance
maxround=20;  %max iteration round
B=X;
comper=1;
for round = 1:maxround
   oldcom=comper;
   for k=1:3
       B=X;
        %% compute the iteration TTM
        for j=1:3
            if(j~=k)
                unfoldA=double(A(:,:,j));
                B=ttm(B,unfoldA',j);
            end
        end
        %% compute leading orthogonal vector 
        
        %reference
        unfoldB=double(tenmat(B,k));
        [rU,~,~]=svd(unfoldB);
        
        unfoldU=U{k};
        %first guess the U using U got by next iteration
        unfoldB=unfoldU'*(double(tenmat(B,k)));

        %nipals algorithm
        %M=unfoldB*unfoldB';
        %[~,newU]=nipals(M,order);
        
        %sys jacobi algorithm
        if(sort==1)
            
        [newU,~,~,sysround]=sysjac(unfoldB,unfoldU',1);
        fprintf(' k= %2d, sysround=%3d\n',k,sysround);
        
        %ring jacobi algorithm
        [newU,~,~,sysround]=ringjacob(unfoldB,unfoldU',1);
        fprintf(' k= %2d, ringround=%3d\n',k,sysround);
            
        A(:,:,k)=newU(:,1:32);
        else
            %% direct iteration way
           A(:,:,k)=ite_svd(double(tenmat(B,k)),order,15);
        end
        %errorU=(norm(A(:,:,k))-norm(rU(:,1:32)))/norm(rU(:,1:32))

        U{k}=newU;

        
   end
   %% compute the error
    core=ttm(B,A(:,:,k),k,'t');
   %{
    for i=1:3
        unfoldA=double(A(:,:,i));
        core=ttm(core,unfoldA,i);
    end
    %}
   %error computation
      error=sqrt(normX^2-norm(core)^2);
      comper=error/normX;
      comdelta=-(oldcom-comper);
      if(sort==1)
      fprintf(' overall_round %2d: comp_error = %e comp_delta = %7.1e\n', round, comper, comdelta);
      end
      if(abs(comdelta)<=comptol)
         break
      end
end
    %% output
        OU=cell(3,1);
    for i=1:3
        OU{i}=double(A(:,:,i));
    end
    T=ttensor(core,OU);
end
