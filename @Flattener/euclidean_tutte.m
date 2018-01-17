     function euclidean_tutte(obj)
        tid=tic;
            cons=PosConstraints(length(obj.M_cut.V));

            obj.flat_V(obj.cut_cone_inds,:)=obj.P;
            
           
          
            obj.flat_V(obj.cut_cone_inds,:)=obj.solver.P;
            
            for j=1:length(obj.M_cut.pathPairs)
                %take the two paths
                path1=obj.M_cut.pathPairs{j}(:,1);
                path2=obj.M_cut.pathPairs{j}(:,2);
                
                Mcur=obj.M{j};
                %indices of previous and next point
                
                %if this is the first iteration, set the position of the
                %boundary to a triangle
                
                x1=linspace(obj.flat_V(path1(1),1),obj.flat_V(path1(end),1),length(path1));
                y1=linspace(obj.flat_V(path1(1),2),obj.flat_V(path1(end),2),length(path1));
                obj.flat_V(path1,:)=[x1' y1'];
                for t=1:length(path1)
                    cons.addConstraint(path1(t),1,[x1(t); y1(t)]);
                end
                x2=linspace(obj.flat_V(path2(1),1),obj.flat_V(path2(end),1),length(path2));
                y2=linspace(obj.flat_V(path2(1),2),obj.flat_V(path2(end),2),length(path2));
                obj.flat_V(path2,:)=[x2' y2'];
                for t=1:length(path2)
                    cons.addConstraint(path2(t),1,[x2(t); y2(t)]);
                end
            end
            %             end
            
            
            
            
            %                 L = compute_combinatorial_laplacian( triangulation2adjacency(TR.ConnectivityList) );
            
            
            
            %                 L=mean_value_laplacian(obj.M_cut.V,obj.M_cut.T);
            %             L=compute_mvc_laplacian(obj.M_cut.V,obj.M_cut.T);
            RealL=sparse(size(obj.L,1)*2,size(obj.L,2)*2);
            RealL(1:2:end,1:2:end)=obj.L;
            RealL(2:2:end,2:2:end)=obj.L;
            L=RealL;
            
            
            %                 fprintf('compute: %f seconds\n',toc(tid));
            
            tidc=tic;
            x=computeFlattening(cons.A,cons.b,L);
            %                 fprintf('lin solve: %f seconds\n',toc(tidc));
            tid=tic;
            
            X=x(1:2:end);
            Y=x(2:2:end);
            obj.flat_V=[X Y];
            obj.times.euclidean_tutte=toc(tid);        
end