function expandedX=expandX(obj,X)
X=X/obj.radius;
            onedim=size(X,2)==1;
            if onedim
                X=[X(1:2:end) X(2:2:end)];
            end
            expandedX=obj.default_X;%zeros(obj.full_len_X,2);
            expandedX(obj.freeVertices,:)=X;
            expandedX=obj.correctX(expandedX);
            if ~onedim
                expandedX=[expandedX(1:2:end) expandedX(2:2:end)];
            end
        end