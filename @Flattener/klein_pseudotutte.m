function klein_pseudotutte(obj)
    tid=tic;
    if isempty(obj.cut_cone_inds)
        obj.init();
    end
    V=obj.hyperbolic_convex_boundary();
    cons=PosConstraints(length(obj.M_cut.V));
    inds=find(~isnan(V(:,1)));
    for t=1:length(inds)
        cons.addConstraint(inds(t),1,[ V(inds(t),1) V(inds(t),2)]);
    end


    RealL=sparse(size(obj.L,1)*2,size(obj.L,2)*2);
    RealL(1:2:end,1:2:end)=obj.L;
    RealL(2:2:end,2:2:end)=obj.L;
    L=RealL;



    x=computeFlattening(cons.A,cons.b,L);

    X=x(1:2:end);
    Y=x(2:2:end);
    V=[X Y];
    obj.flat_V=klein_to_poincare(V);

    obj.times.klein_pseudotutte=toc(tid);
end