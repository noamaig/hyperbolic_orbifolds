  function X=correctX(obj,X)
  
            if size(X,2)==1
                X=[X(1:2:end) X(2:2:end)];
            end
            
            X(~obj.freeVertices,:)=obj.default_X(~obj.freeVertices,:);
            for i=1:length(obj.pathPairs)
                p=obj.pathPairs{i};
                X(p(:,2),:)=obj.mobius_trans{i}.map(X(p(:,1),:));
            end
            X(obj.cone_inds,:)=obj.P;
            
            X=X';
            X=X(:);
        end
      