function reducedX=reduceX(obj,X)
X=X*obj.radius;
            onedim=size(X,2)==1;
            if onedim
                reducedX=[X(1:2:end) X(2:2:end)];
            else
                reducedX=X;
            end
            reducedX=reducedX(obj.freeVertices,:);
            if onedim
                reducedX=reducedX';
                reducedX=reducedX(:);
            end
        end
        