function vprops=face2vertex_int(obj,property,areas)
            inds=reshape(obj.M_orig.F,1,3*length(obj.M_orig.F))';
            areas=kron(areas,[1; 1; 1]);
            props=kron(property,[1; 1; 1]);
            vareas=accumarray(inds,areas);
            vprops=accumarray(inds,areas.*props);
            vprops=vprops./vareas;
        end