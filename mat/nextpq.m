function [P,Q] = nextpq(P,Q,i,pol)
%NEXTPQ decide next P,Q in ring jacobi array
    a=2*size(P,2);
    if(pol==1) %forward sweep
        sn=ceil((i)/2);
        temp=Q(ceil(i/2));
        IQ=[Q(1:sn),P(sn),Q(sn+1:a/2)];
        IQ=[IQ(a/2+1),IQ(1:a/2)];
        if(sn>=a/2)
            Q=IQ(1:sn);
        else
        Q=[IQ(1:sn),IQ(sn+2:a/2+1)];
        end
        IP=[P(1:sn),temp,P(sn+1:a/2)];
        if(sn<2)
        P=IP(2:a/2+1);
        else
        P=[IP(1:sn-1),IP(sn+1:a/2+1)];
        end
        
    else %backward sweep
        sn=a/2-ceil((i-2)/2);
        temp=P(sn);
        IP=[P(1:sn-1),Q(sn),P(sn:a/2)];
        IP=[IP(2:a/2+1),IP(1),];
        if(sn==1)
            P=IP(2:a/2+1);
        else
        P=[IP(1:sn-1),IP(sn+1:a/2+1)];
        end
        
        IQ=[Q(1:sn-1),temp,Q(sn:a/2)];
        if(sn==a/2)
        Q=IQ(1:a/2);
        else
        Q=[IQ(1:sn),IQ(sn+2:a/2+1)];
        end
    end
    
end

