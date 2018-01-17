function X=gd(obj,X)
            lastX=inf(size(X));
            lastO=inf;
            iter=0;
            while(true)     
                iter=iter+1;
                if iter>200
                    break;
                end
                X=obj.reduceX(X);
                [X,O]=obj.gd_step(X);
                X=obj.expandX(X);
                X=obj.correctX(X);
                
                Xchange=norm(X-lastX);
%                 if mod(iter,50)==1
                 
                    fprintf('change in X: %g, objective: %g\n',Xchange,lastO-O);
                 
%                 end
                if Xchange<1e-16
                    break;
                end
                
                lastX=X;
                lastO=O;
            end
            
        end
        