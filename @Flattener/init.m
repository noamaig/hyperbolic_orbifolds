function init(obj)
tid=tic;
TR=triangulation(obj.M_orig.F',obj.M_orig.V');
bdry=TR.freeBoundary();
if ~isempty(bdry)
    bdry=bdry(:,1);    
else
    if isempty(obj.M_cut)
        
        obj.cut();
    end
    startP=obj.M_cut.Old2New{obj.uncut_cone_inds(1)};
    assert(length(startP)==1);
    TR=triangulation(obj.M_cut.T,obj.M_cut.V);
    pathEnds=[];
    for i=1:length(obj.M_cut.pathPairs)
        pathEnds=[pathEnds obj.M_cut.pathPairs{i}([1 end],:)];
    end
    pathEnds=unique(pathEnds);
    
    all_binds = TR.freeBoundary();
    assert(all(all_binds(:,2)==all_binds([2:end,1],1)));
    all_binds=all_binds(:,1);
    ind=find(all_binds==startP);
    all_binds=all_binds([ind:end,1:ind-1]);
    obj.cut_cone_inds=find(ismember(all_binds,pathEnds));
    obj.cut_cone_inds=all_binds(obj.cut_cone_inds);
    obj.cut_cone_inds=obj.cut_cone_inds(end:-1:1);
    obj.cut_cone_inds=[obj.cut_cone_inds(end); obj.cut_cone_inds(1:end-1)];
    
    assert(length(obj.cut_cone_inds)==length(obj.P));
end

L = cotmatrix(obj.M_cut.V,obj.M_cut.T);

obj.L=L;

clamp=1e-2;
inds=find(obj.L~=0);
temp=obj.L(inds)<clamp;
inds=inds(temp);
obj.L(inds)=clamp;
obj.L(sub2ind(size(obj.L),1:length(obj.L),1:length(obj.L)))=0;
obj.L(sub2ind(size(obj.L),1:length(obj.L),1:length(obj.L)))=-sum(obj.L);
obj.times.init=toc(tid);
end