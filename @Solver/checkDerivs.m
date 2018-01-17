  function check_derivs(obj,X)
            %transforming X for easiness
            orgX=X;
            for iter=1:20
                X=orgX;
                X=obj.correctX(X);
                X=[X(1:2:end) X(2:2:end)];
              
                step_size=1e-8;
                delta=(rand(size(X))-0.5)*step_size;
                %the offset vector
                Xd=X+delta;
                %making sure the offset vector complies to the given
                %constraints:
                %1. mobius
                for i=1:length(obj.pathPairs)
                    p=obj.pathPairs{i};
                    Xd(p(:,2),:)=obj.mobius_trans{i}.map(Xd(p(:,1),:));
                end
                %2. cones
                Xd(obj.cones,:)=X(obj.cones,:);
                %             Xd(obj.inner_vertices,:)=X(obj.inner_vertices,:);
                Xd=Xd';
                Xd=Xd(:);
                X=X';
                X=X(:);
                delta=Xd-X;
                
                o=obj.objectives(X);
                
                [o,g]=obj.objective_and_grad(X);
                od=obj.objective_and_grad(X+delta);
                n_deriv=sum(g.*delta);
                fd_deriv=(sum(od.^2)-sum(o.^2));
                d(iter)=(n_deriv/fd_deriv)
            end
        end