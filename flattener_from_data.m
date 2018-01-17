function flattener = flattener_from_data( fname )
%load saved file to a flattener object
s=load(fname);
if isfield(s,'flattener')
    s=s.flattener;
    V=s.M_orig.V;
    T=s.M_orig.F;
    inds=s.uncut_cone_inds;
else
    [V,T]=read_off(['models/' s.mesh_name]);
    inds=s.inds;
end
if all(V(3,:)==0)
    V(3,:)=[];
end
V=V/mean(std(V'));
V=bsxfun(@minus,V',mean(V'))';
flattener=Flattener(V,T,inds,'isdisc',size(V,1)==2);
flattener.mesh_name=s.mesh_name;

end

