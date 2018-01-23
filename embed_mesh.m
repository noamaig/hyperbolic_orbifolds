function [flattener] = embed_mesh(V,T,inds,isdisc,forceOrder)
%main function for mesh embedding. 
%V - vertices
%T - triangles
%inds - the indices of the vertices that are cones
%isdisc - boolean flag, true iff the mesh is a disk (otherwise must be sphere).
%force order - if true use given cone order. if not given or false, use a
%Travelling Salesmen approximation to compute the ordering which will
%generate the shortest cut graph.
if nargin<5
    forceOrder=false;
end
flattener=Flattener(V,T,inds,'isdisc',isdisc);
flattener.mesh_name='mesh';
if isdisc
    flattener.disk_orbifold();
else
    if ~forceOrder
        flattener.orderTS();
    end
    flattener.flatten_orbifold(true);
end
end

