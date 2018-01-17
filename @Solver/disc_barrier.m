 function [f,g]=disc_barrier(obj,X)
            barrier_fac=1e+8;
            X=[X(1:2:end) X(2:2:end)];
            assert(size(X,2)==2);
            n=sum(X.^2,2);
            n=max(n-1,0);
            f=barrier_fac*sum(exp(n.^2)-1);
            g=zeros(size(X));
            %d(n)=d(e^(x^2+y^2-1))=e^(x^2+y^2-1)2x,
            g=barrier_fac*2*[n.*X(:,1) n.*X(:,2)];
            g=g';g=g(:);
        end