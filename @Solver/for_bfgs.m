function [o,g]=for_bfgs(obj,X)
orgX=X;
X=obj.expandX(X);
if size(X,2)==1
    X=[X(1:2:end) X(2:2:end)];
end
for i=1:length(obj.pathPairs)
    p=obj.pathPairs{i};
    X(p(:,2),:)=obj.mobius_trans{i}.map(X(p(:,1),:));
end
X(obj.cone_inds,:)=obj.P;
d=sum(X.^2,2);
%             if any(d>1)
%                 o=inf;
%                 g=zeros(length(X)*2,1);
%                 return;
%             end
X=X';
X=X(:);
%             obj.check_derivs(X);
%             o=obj.objectives(X);
%             j=obj.GN_jacobian(X);
%             g=j'*o;
[d,g,w]=obj.objective_and_grad(X);
%next step is heuristic division by the metric for
%"precondtioning"
if 0 %obj.precondition
    mtemp=obj.metric(X);
    m=zeros(length(mtemp)*2,1);
    m(1:2:end)=mtemp;
    m(2:2:end)=mtemp;
    g=g./m;
end
o=sum(w.*(d.^2));
[fb,gb]=obj.disc_barrier(X);
g=g+gb;
g=obj.reduceX(g);
%!!! addition
if length(orgX)<500
    o=obj.karcherPerVertex(X);
    o=obj.reduceX([o o]);
    o=sum(o(:,1));
end
o=o+fb;
end
