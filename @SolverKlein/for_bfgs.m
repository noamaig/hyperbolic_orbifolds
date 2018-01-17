function [o,g]=for_bfgs(obj,X)
% orgX=X;
% X=obj.expandX(X);
% if size(X,2)==1
%     X=[X(1:2:end) X(2:2:end)];
% end
% for i=1:length(obj.pathPairs)
%     p=obj.pathPairs{i};
%     X(p(:,2),:)=obj.mobius_trans{i}.map(X(p(:,1),:));
% end
% X(obj.cone_inds,:)=obj.P;
% d=sum(X.^2,2);
%             if any(d>1)
%                 o=inf;
%                 g=zeros(length(X)*2,1);
%                 return;
%             end
% X=X';
% X=X(:);
%             obj.check_derivs(X);
%             o=obj.objectives(X);
%             j=obj.GN_jacobian(X);
%             g=j'*o;
[hd,g,w]=obj.objective_and_grad(X);
%next step is heuristic division by the metric for
%"precondtioning"
if obj.precondition
%     error('set metric');
    %
    %first term of metric is I/(1-|P|^2)
    %second term is 
    %[X^2 XY; XY Y^2]/(1-|P|^2)^2
    %
    x=X(1:2:end);
    y=X(2:2:end);
    nrm=x.^2+y.^2;
    denom=1-nrm;
    a=1./denom+(x.^2)./denom.^2;
    b=(x.*y)./denom.^2;
    c=b;
    d=1./denom+(y.^2)./denom.^2;
    
    dtr=a.*d-b.*c;
    ai=d./dtr;
    bi=-b./dtr;
    ci=-c./dtr;
    di=a./dtr;
%     warning('trying new heac');
%     ai=denom.^2;%1./sqrt(dtr);
%     bi=denom.^2;%1./sqrt(dtr);
%     ci=denom.^2;%1./sqrt(dtr);
%     di=denom.^2;%1./sqrt(dtr);
    gx=g(1:2:end);
    gy=g(2:2:end);
    
    newgx=ai.*gx+bi.*gy;
    newgy=ci.*gx+di.*gy;
    
    g=[newgx newgy];
    g=g';
    g=g(:);
    
    
end
o=sum(w.*(hd.^2));
[fb,gb]=obj.disc_barrier(X);
g=g+gb;
% g=obj.reduceX(g);

o=o+fb;
end
