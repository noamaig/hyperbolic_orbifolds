function [vcol,ncol]=computeVCol(obj,uncut,varargin)
if uncut
    vcol=obj.V;
else
    vcol=obj.M_cut.V;
end
if size(vcol,2)==2
    vcol=[vcol ones(length(vcol),1)];
end
vcol=vcol-repmat(min(vcol),length(vcol),1);
vcol=vcol./repmat(max(vcol),length(vcol),1);
end