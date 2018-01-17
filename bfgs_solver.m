classdef bfgs_solver <handle
    properties
        inds;
        W;
        Wright;
        Winterior;
        WRightToRight;
        I;
        J;
        Iinterior;
        Jinterior;
        inverse_mobius_a;
        inverse_mobius_b;
        inverse_mobius_c;
        inverse_mobius_d;
        Iright;
        Jright;
        P;
        cone_inds;
        pathPairs={};
        mobius_trans={};
        polygon_inds={};
        polygon_P={};
        inner_vertices;
        left_side_vertices;
        right_side_vertices;
        right_to_left;
        cones;
        adj;
        path_ind;
        IRightToRight;
        JRightToRight;
        old_p;
        cones_cons;
        L;
        fix_gradient_tol;
    end
    methods
        
        function H=hessFun(obj,X,lambda)
            H=obj.L;
        end
        function [f,g]=disc_barrier(obj,X)
            barrier_fac=1e+8;
            X=[X(1:2:end) X(2:2:end)];
            assert(size(X,2)==2);
            n=sum(X.^2,2);
            n=max(n-1,0);
            f=barrier_fac*sum(exp(n.^2)-1);
            g=zeros(size(X));
            %d(n)=d(e^(x^2+y^2-1))=e^(x^2+y^2-1)2x,
            g=barrier_fac*2*[n.*X(:,1) n.*X(:,2)];
            g=g';g=g(:);
        end
        function x=solve_bfgs(obj,X)
            if obj.fix_gradient_tol>0
                [~,g]=obj.objective_and_grad(X);
                g=reshape(g,2,round(length(g)/2))';
                x=reshape(X,2,round(length(X)/2))';
                n=sqrt(sum(g.^2,2));
                inds=find(n<obj.fix_gradient_tol);
                inds=setdiff(inds,[find(obj.right_side_vertices);obj.cone_inds]);
                            obj.cones_cons=PosConstraints(length(obj.adj));

                for i=1:length(inds)
               
                obj.cones_cons.addConstraint(inds(i),1,x(i,:));
                
            
                end
                
            end
            options =optimoptions('fmincon');
            options.MaxFunEvals=inf;
            options.Display = 'iter-detailed';
            options.Diagnostics = 'on';
            %             params{SOLVER_NO}.options.DerivativeCheck = 'on';
            
            options.GradObj = 'on';
            options.MaxIter = 50000;
            options.TolFun=0;1e-10;
            options.TolX=1e-32;
            if isempty(obj.L)
                options.Hessian = {'lbfgs',1};
            else
                options.Hessian = 'user-supplied';
                options.HessFcn=@(X,lambda)obj.hessFun(X,lambda);
            end
%             assert(norm(obj.cones_cons.A*X-obj.cones_cons.b)<1e-10);
            x = fmincon(@(X)obj.for_bfgs(X),X,[],[],obj.cones_cons.A,obj.cones_cons.b,[],[],[],options);
%             options=[];
%             options.optTol=1e-10;
%             options.progTol=1e-12;
%             options.MaxFunEvals=1e+8;
%             options.MaxIter=500;
%             x = minFunc(@(X)obj.for_bfgs(X),X,options);
            x=obj.correctX(x);
        end
        
        
        function m=metric(obj,X)
            m=4./(1-X(1:2:end).^2-X(2:2:end).^2).^2;
        end
        
        
        function [o,g]=for_bfgs(obj,X)
            if size(X,2)==1
                X=[X(1:2:end) X(2:2:end)];
            end
%             for i=1:length(obj.pathPairs)
%                 p=obj.pathPairs{i};
%                 X(p(:,2),:)=obj.mobius_trans{i}.map(X(p(:,1),:));
%             end
            X(obj.cone_inds,:)=obj.P;
            d=sum(X.^2,2);
%             if any(d>1)
%                 o=inf;
%                 g=zeros(length(X)*2,1);
%                 return;
%             end
            X=X';
            X=X(:);
            %             obj.check_derivs(X);
            %             o=obj.objectives(X);
            %             j=obj.GN_jacobian(X);
            %             g=j'*o;
            [d,g,w]=obj.objective_and_grad(X);
            %next step is heuristic division by the metric for
            %"precondtioning"
            mtemp=obj.metric(X);
            m=zeros(length(mtemp)*2,1);
            m(1:2:end)=mtemp;
            m(2:2:end)=mtemp;
            g=g./m;
            
            o=sum(w.*(d.^2));
            [fb,gb]=obj.disc_barrier(X);
            g=g+gb;
            o=o+fb;
        end
        
        function X=correctX(obj,X)
            if size(X,2)==1
                X=[X(1:2:end) X(2:2:end)];
            end
            for i=1:length(obj.pathPairs)
                p=obj.pathPairs{i};
                X(p(:,2),:)=obj.mobius_trans{i}.map(X(p(:,1),:));
            end
            X(obj.cone_inds,:)=obj.P;
            
            X=X';
            X=X(:);
        end
        function check_derivs(obj,X)
            %transforming X for easiness
            orgX=X;
            for iter=1:20
                X=orgX;
                X=obj.correctX(X);
                X=[X(1:2:end) X(2:2:end)];
                
                %             hold on;
                %             for i=1:length(obj.IRightToRight)
                %                 line(X([obj.IRightToRight(i);obj.JRightToRight(i)],1),X([obj.IRightToRight(i);obj.JRightToRight(i)],2));
                %             end
                %             pause
                %an offset
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
                %             assert(max(abs(d))<1e-6);
            end
        end
        function obj=bfgs_solver(adj,pathPairs,mobius_trans,cone_inds,P,varargin)
            parser = inputParser;
            parser.addOptional('fix_gradient_tol',0,@isnumeric);
            parser.parse(varargin{:});
            obj.P=P;
            obj.adj=adj~=0;
            [obj.I,obj.J]=find(adj);
            od=obj.I~=obj.J;
            obj.I=obj.I(od);
            obj.J=obj.J(od);
            W=full(adj(sub2ind(size(adj),obj.I,obj.J)));
            clamp=1e-1;
            W(W<0)=1;
            obj.W=W;
            
            obj.cone_inds=cone_inds;
            
            obj.pathPairs=pathPairs;
            obj.mobius_trans=mobius_trans;
            
            for i=1:length(mobius_trans)
                inverse_mobius=mobius_trans{i}.inverse();
                obj.inverse_mobius_a(i)=inverse_mobius.a;
                obj.inverse_mobius_b(i)=inverse_mobius.b;
                obj.inverse_mobius_c(i)=inverse_mobius.c;
                obj.inverse_mobius_d(i)=inverse_mobius.d;
            end
            obj.inverse_mobius_a=obj.inverse_mobius_a(:);
            obj.inverse_mobius_b=obj.inverse_mobius_b(:);
            obj.inverse_mobius_c=obj.inverse_mobius_c(:);
            obj.inverse_mobius_d=obj.inverse_mobius_d(:);
            obj.inner_vertices=true(length(adj),1);
            obj.left_side_vertices=false(length(adj),1);
            obj.right_side_vertices=false(length(adj),1);
            obj.cones=false(length(adj),1);
            obj.path_ind=-ones(length(adj),1);
            for i=1:length(pathPairs)
                obj.path_ind(obj.pathPairs{i}(:))=i;
                obj.inner_vertices(obj.pathPairs{i}(:))=false;
                obj.left_side_vertices(obj.pathPairs{i}(2:end-1,1))=true;
                obj.right_side_vertices(obj.pathPairs{i}(2:end-1,2))=true;
                obj.cones(obj.pathPairs{i}([1 end],:))=true;
                obj.right_to_left(obj.pathPairs{i}(:,2))=obj.pathPairs{i}(:,1);
            end
            obj.right_to_left=obj.right_to_left';
            
            
            I=obj.I;
            J=obj.J;
            
            %find edges that connect two right sides of seams but do not
            %belong to seam path
            keep=(obj.path_ind(I)>obj.path_ind(J))& ...
                obj.path_ind(J)~=-1 ...
                &obj.right_side_vertices(I)&obj.right_side_vertices(J) ...
                & ~obj.cones(I) & ~obj.cones(J);
            obj.IRightToRight=I(keep);
            obj.JRightToRight=J(keep);
            W=obj.W;
            obj.WRightToRight=W(keep);
            %             assert(all(obj.right_side_vertices(obj.IRightToRight)));
            %             assert(all(obj.right_side_vertices(obj.JRightToRight)));
            %             %discard of edges that are both on right side
            keep=(~obj.right_side_vertices(I))|(~obj.right_side_vertices(J));
            I=I(keep);
            J=J(keep);
            W=W(keep);
            
            %keep only edges that don't end on right side
            %             error('fix for cones');
            keep=~obj.right_side_vertices(J);
            I=I(keep);
            J=J(keep);
            W=W(keep);
            %remove duplicate edges
            newinds=[I J];
            [~,uinds]=unique(sort(newinds,2),'rows');
            I=I(uinds);
            J=J(uinds);
            W=W(uinds);
            obj.I=I;
            obj.J=J;
            obj.W=W;
            %divide to two groups: not right vertices, for which we
            %derive by the first vertex of the edge, and right side
            %vertices for which we transform to left side the interior neighbours and derive by
            %first vertex of the edge. should be what I did? vcheck where
            %fails
            %             assert(all((obj.right_side_vertices(I)~=obj.left_side_vertices(J))|obj.cones(obj.I)));
            right_boundary_edges=obj.right_side_vertices(I);
            %remove all edges emenating from right boundary or cones
            keep=~right_boundary_edges;
            obj.Iinterior=obj.I(keep);
            obj.Jinterior=obj.J(keep);
            obj.Winterior=obj.W(keep);
            %deal with edges emenating from right boundary which are not
            %cones
            keep=right_boundary_edges;
            
            obj.Iright=obj.I(keep);
            obj.Jright=obj.J(keep);
            obj.Wright=obj.W(keep);
            obj.cones_cons=PosConstraints(length(adj));
            for i=1:length(obj.cone_inds)
                obj.cones_cons.addConstraint(obj.cone_inds(i),1,obj.P(i,:));
                
            end
            inds=setdiff(find(obj.right_side_vertices),obj.cone_inds);
            for i=1:length(inds)
                obj.cones_cons.addConstraint(inds(i),1,[0 0]);
            end
            obj.fix_gradient_tol= parser.Results.fix_gradient_tol;
                
            %             assert(length(obj.Iright)+length(obj.Iinterior)==length(obj.I));
        end
        
        
        function hd = karcher_objective(obj, I,J,X )
            assert(size(X,2)==2);
            %making sure mobius trans is correct
            
            
            a=X(I,1);
            b=X(I,2);
            c=X(J,1);
            d=X(J,2);
            delta=2*((a-c).^2+(b-d).^2)./((1-a.^2-b.^2).*(1-c.^2-d.^2));
            hd=acosh(1+delta);
            
            if any(~(isreal(hd)))
                warning('X out of unit disk');
                hd(~isreal(hd))=1;
            end
            
        end
        function o=objectives(obj,X)
            if size(X,2)==1
                X=[X(1:2:end) X(2:2:end)];
            end
            
            interior=obj.karcher_objective(obj.Iinterior,obj.Jinterior,X);
            right=obj.karcher_objective(obj.Iright,obj.Jright,X);
            right_to_right=obj.karcher_objective(obj.IRightToRight,obj.JRightToRight,X);
            o=[interior;right;right_to_right];
        end
        
        function [grada2,gradb2,gradc2,gradd2,hd2]=karcher_grad(obj,a,b,c,d,w)
            
            
        f=2*(a.^2+b.^2+c.^2+d.^2-2*a.*c-2*b.*d);
        g=1-a.^2-b.^2-c.^2-d.^2+a.^2.*c.^2+a.^2.*d.^2+b.^2.*c.^2+b.^2.*d.^2;
        delta2=f./g;
        hd2=real((acosh(1+delta2)));
        f_a=4*a-4*c;
        f_b=4*b-4*d;
        f_c=4*c-4*a;
        f_d=4*d-4*b;
        g_a=2*a.*(c.^2+d.^2-1);
        g_b=2*b.*(c.^2+d.^2-1);
        g_c=2*c.*(a.^2+b.^2-1);
        g_d=2*d.*(a.^2+b.^2-1);
        
        grad_delta_a=(g.*f_a-f.*g_a)./(g.^2);
        grad_delta_b=(g.*f_b-f.*g_b)./(g.^2);
        grad_delta_c=(g.*f_c-f.*g_c)./(g.^2);
        grad_delta_d=(g.*f_d-f.*g_d)./(g.^2);
        
        fac=(sqrt(delta2.^2+2*delta2));
        grada2=2*w.*hd2.*grad_delta_a./fac;
        gradb2=2*w.*hd2.*grad_delta_b./fac;
        gradc2=2*w.*hd2.*grad_delta_c./fac;
        gradd2=2*w.*hd2.*grad_delta_d./fac;
        grada2(hd2==0)=0;
        gradb2(hd2==0)=0;
        gradc2(hd2==0)=0;
        gradd2(hd2==0)=0;
        
%         ab=(1-a.^2-b.^2).^2;
%         cd=(1-c.^2-d.^2).^2;
%         
%         grada2=grada2.*ab;
%         gradb2=gradb2.*ab;
%         
%         gradc2=gradc2.*cd;
%         gradd2=gradd2.*cd;
        
%         hd2=real(hd2);

            
            
            

        end
        
        
        function [O,G,w]=objective_and_grad(obj,X)
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
        
        
        
        
    end
end

