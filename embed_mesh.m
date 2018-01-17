function [flattener] = embed_mesh(V,T,inds,isdisc)
flattener=Flattener(V,T,inds,'isdisc',isdisc);
flattener.mesh_name='mesh';
if isdisc
flattener.disk_orbifold();
else
    flattener.flatten_orbifold(true);
end
% flattener.computeDistortion(true);
end

