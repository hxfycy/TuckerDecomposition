function [S,U,error] = jacobiunit(pl,ql,S,U,sort)
%JACOBIUNIT One single jacobi unit 
n=size(S,1);
%do the rotation
        tmpb=S(pl,:);
        tmpu=U(pl,:);
    gama=S(pl,:)*S(ql,:)';
    alpha=S(pl,:)*S(pl,:)';
    beta=S(ql,:)*S(ql,:)';
    
        % change two column's term¡ª¡ªThis work can be implemented in
        % hardware realy easily
        if(alpha<beta&&pl<ql&&sort==1)
            S(pl,:)=S(ql,:);
            S(ql,:)=tmpb;
            U(pl,:)=U(ql,:);
            U(ql,:)=tmpu;
            tmpb=S(pl,:);
            tmpu=U(pl,:);
            beta=S(ql,:)*S(ql,:)';
            alpha=S(pl,:)*S(pl,:)';
        end
    sita=atan(2*gama/(beta-alpha));
    sita=sita/2;
    error=abs(gama)/sqrt(alpha*beta);
    S(pl,:)=tmpb*cos(sita)-S(ql,:)*sin(sita);
    S(ql,:)=tmpb*sin(sita)+S(ql,:)*cos(sita);
    U(pl,:)=cos(sita)*tmpu-sin(sita)*U(ql,:);
    U(ql,:)=tmpu*sin(sita)+U(ql,:)*cos(sita);
       
end

