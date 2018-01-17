    
        function [O,G,w]=objective_and_grad(obj,X)
            assert(size(X,2)==1);
            assert(~any(isnan(X(:))));
            %given a vector s, return the jacboain G of the karcher
            %objectives, that is
            %G_ij is the hyperbolic distance between the vertices of the
            %i'th edge, derived according to the j'th variable
            %karcher
            
            
            I=obj.Iinterior;
            J=obj.Jinterior;
            a=X(I*2-1);
            b=X(I*2);
            c=X(J*2-1);
            d=X(J*2);
            [grada,gradb,gradc,gradd,O]=obj.karcher_grad(a,b,c,d,obj.Winterior);
            w=obj.Winterior;
            
            Ja=I*2-1;
            Jb=I*2;
            Jc=J*2-1;
            Jd=J*2;
            
            J_k=[Ja;Jb;Jc;Jd];
            
            V_k=[grada;gradb;gradc;gradd];
            
            %now taking only edges which are from right boundary to
            %interior and not cones
            I=obj.Iright;
            J=obj.Jright;
            
            %go over all paths and build the 1 ring of each left vertex and
            %compute gradient w.r.t. it.
            
            %the index of the vertex on right side of the cut
            right_ind_I=I;
            %                 assert(obj.right_side_vertices(right_ind_I));
            %the left counter point of the right vertex
            left_ind_I=obj.right_to_left(right_ind_I);
            %the interior vertex connected to vertex on right side of
            %cut
            right_ind_J=J;
            %                 assert(obj.r(right_ind_J));
            %the path the right vertex is on
            path_ind=obj.path_ind(right_ind_I);
            %mobius trans taking left side to right side
            
            %take the interior vertex using the mobius trans to
            %corrsponding place for right side
            tomap=X(right_ind_J*2-1)+1i*X(right_ind_J*2);
            ma=obj.inverse_mobius_a(path_ind);
            mb=obj.inverse_mobius_b(path_ind);
            mc=obj.inverse_mobius_c(path_ind);
            md=obj.inverse_mobius_d(path_ind);
            Y=(ma.*tomap+mb)./(mc.*tomap+md);
            Y=[real(Y) imag(Y)];
            %                 assert(norm(M.map(Y)-[X(right_ind_J*2-1) X(right_ind_J*2)])<1e-10);
            %now compute gradient.
            Ja=left_ind_I*2-1;
            Jb=left_ind_I*2;
            Jc=right_ind_J*2-1;
            Jd=right_ind_J*2;
            %derivatice w.r.t the 1st vertex of edge which is on the
            %right side
            a=X(Ja);
            b=X(Jb);
            c=Y(:,1);
            d=Y(:,2);
            [grada,gradb,~,~,o]=obj.karcher_grad(a,b,c,d,obj.Wright);
            w=[w;obj.Wright];
            O=[O;o];
            %derivatice w.r.t the 2nd vertex of edge which is *NOT* on the
            %right side
            a=X(right_ind_I*2-1);
            b=X(right_ind_I*2);
            c=X(right_ind_J*2-1);
            d=X(right_ind_J*2);
            [~,~,gradc,gradd]=obj.karcher_grad(a,b,c,d,obj.Wright);
            %             O=[O;o];
            %add gradient to vector
            V_k=[V_k;grada;gradb;gradc;gradd];
            %indices of relevant entries - we imagine the tranformed
            %inner vertex is connected to the LEFT side (since we want
            %to ignore right side completely during optimization)
            
            J_k=[J_k;Ja;Jb;Jc;Jd];
            if ~isempty(obj.IRightToRight)
                I=obj.IRightToRight;
                J=obj.JRightToRight;
                
                %go over all paths and build the 1 ring of each left vertex and
                %compute gradient w.r.t. it.
                
                %the index of the vertex on right side of the cut
                right_ind_I=I;
                right_ind_J=J;
                %the left counter point of the right vertex
                left_ind_I=obj.right_to_left(right_ind_I);
                left_ind_J=obj.right_to_left(right_ind_J);
                %the interior vertex connected to vertex on right side of
                %cut
                
                %the path the right vertex is on
                path_ind_I=obj.path_ind(right_ind_I);
                path_ind_J=obj.path_ind(right_ind_J);
                %mobius trans taking left side to right side
                
                tomap_I=X(right_ind_J*2-1)+1i*X(right_ind_J*2);
                ma_I=obj.inverse_mobius_a(path_ind_I);
                mb_I=obj.inverse_mobius_b(path_ind_I);
                mc_I=obj.inverse_mobius_c(path_ind_I);
                md_I=obj.inverse_mobius_d(path_ind_I);
                Y_I=(ma_I.*tomap_I+mb_I)./(mc_I.*tomap_I+md_I);
                Y_I=[real(Y_I) imag(Y_I)];
                
                tomap_J=X(right_ind_I*2-1)+1i*X(right_ind_I*2);
                ma_J=obj.inverse_mobius_a(path_ind_J);
                mb_J=obj.inverse_mobius_b(path_ind_J);
                mc_J=obj.inverse_mobius_c(path_ind_J);
                md_J=obj.inverse_mobius_d(path_ind_J);
                Y_J=(ma_J.*tomap_J+mb_J)./(mc_J.*tomap_J+md_J);
                Y_J=[real(Y_J) imag(Y_J)];
                
                
                %take the interior vertex using the mobius trans to
                %corrsponding place for right side
                %                 Y_I=M_J.inverse().map([X(right_ind_I*2-1) X(right_ind_I*2)]);
                %                 Y_J=M_I.inverse().map([X(right_ind_J*2-1) X(right_ind_J*2)]);
                
                %now compute gradient.
                Ja=left_ind_I*2-1;
                Jb=left_ind_I*2;
                Jc=left_ind_J*2-1;
                Jd=left_ind_J*2;
                %derivatice w.r.t the 1st vertex of edge which is on the
                %right side
                a=X(Ja);
                b=X(Jb);
                c=Y_I(:,1);
                d=Y_I(:,2);
                [grada,gradb,~,~,o]=obj.karcher_grad(a,b,c,d,obj.WRightToRight);
                w=[w;obj.WRightToRight];
                O=[O;o];
                %derivatice w.r.t the 2nd vertex of edge which is *NOT* on the
                %right side
                a=Y_J(:,1);
                b=Y_J(:,2);
                c=X(Jc);
                d=X(Jd);
                [~,~,gradc,gradd]=obj.karcher_grad(a,b,c,d,obj.WRightToRight);
                %                 O=[O;o];
                %add gradient to vector
                V_k=[V_k;grada;gradb;gradc;gradd];
                %indices of relevant entries - we imagine the tranformed
                %inner vertex is connected to the LEFT side (since we want
                %to ignore right side completely during optimization)
                
                J_k=[J_k;Ja;Jb;Jc;Jd];
                
                %             end
            end
            if max(J_k)<length(obj.adj)
                J_k=[J_k;length(obj.adj)];
                V_k=[V_k;0];
            end
            
            G = accumarray(J_k, V_k);
            inds=find(obj.cones);
            G(inds*2-1)=0;
            G(inds*2)=0;
            %make sure the "redundant" vertices (cones + right side) have
            %no effect on jacobian.
            %             rs=find(obj.right_side_vertices);
            %             cones=find(obj.cones);
            %             assert(all(all(G(rs*2-1)==0)));
            %             assert(all(all(G(rs*2)==0)));
            %             assert(all(all(G(cones*2-1)==0)));
            %             assert(all(all(G(cones*2)==0)));
        end
    