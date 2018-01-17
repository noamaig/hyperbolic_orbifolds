function cut(obj)
    tid=tic;
    root=1;
    fixedPairs=[1:length(obj.uncut_cone_inds)-1;2:length(obj.uncut_cone_inds)]';
    tree=sparse(fixedPairs(:,1),fixedPairs(:,2),1,length(obj.uncut_cone_inds),length(obj.uncut_cone_inds));
    cutter=TreeCutter(obj.M_orig.V',obj.M_orig.F',tree,obj.uncut_cone_inds,root);%+(rand(size(obj.M_orig.V'))*20)*diag([1 100 100])+rand(size(obj.M_orig.V'))*20
    cutter.cutTree();
    obj.setCutMesh(cutter);
    obj.times.cut=toc(tid);
end
