    
        function [O,G,w]=objective_and_grad(obj,X)
            assert(size(X,2)==1);
            assert(~any(isnan(X(:))));
            %given a vector s, return the jacboain G of the karcher
            %objectives, that is
            %G_ij is the hyperbolic distance between the vertices of the
            %i'th edge, derived according to the j'th variable
            %karcher
            
            
            I=obj.I;
            J=obj.J;
            a=X(I*2-1);
            b=X(I*2);
            c=X(J*2-1);
            d=X(J*2);
            [grada,gradb,gradc,gradd,O]=obj.karcher_grad(a,b,c,d,obj.W);
            w=obj.W;
            
            Ja=I*2-1;
            Jb=I*2;
            Jc=J*2-1;
            Jd=J*2;
            
            J_k=[Ja;Jb;Jc;Jd];
            
            V_k=[grada;gradb;gradc;gradd];
            
           
            if max(J_k)<length(obj.adj)
                J_k=[J_k;length(obj.adj)];
                V_k=[V_k;0];
            end
            
            G = accumarray(J_k, V_k);
            assert(~any(isnan(G)));
%             inds=find(obj.cones);
%             G(inds*2-1)=0;
%             G(inds*2)=0;
            %make sure the "redundant" vertices (cones + right side) have
            %no effect on jacobian.
            %             rs=find(obj.right_side_vertices);
            %             cones=find(obj.cones);
            %             assert(all(all(G(rs*2-1)==0)));
            %             assert(all(all(G(rs*2)==0)));
            %             assert(all(all(G(cones*2-1)==0)));
            %             assert(all(all(G(cones*2)==0)));
        end
    