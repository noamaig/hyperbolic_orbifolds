function [inds] = knn(X,Y)
    inds=arrayfun(@(i)nn(X,Y(i,:)),1:size(Y,1));
end
function ind=nn(X,y)
    b=sum(bsxfun(@minus,X,y).^2,2);
    [~,ind]=min(b);
end

