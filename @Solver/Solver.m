classdef Solver <handle
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
        L;
        fix_gradient_tol;
        freeVertices;
        full_len_X;
        default_X;
        Wmat;
        precondition=true;
        radius=1;
        
        % aux for faster computations
        J_k_accumarray_mat;
        V_temp;
    end
    methods
        
        
       
        
        
       
        
        
        
      
        function obj=Solver(adj,pathPairs,mobius_trans,cone_inds,P,varargin)
            parser = inputParser;
            parser.addOptional('fix_gradient_tol',0,@isnumeric);
            parser.addOptional('fix_inds',[],@(X)all(isnumeric(X)));
            parser.addOptional('default_x',zeros(length(adj),2));
            parser.parse(varargin{:});
                         obj.full_len_X=length(adj);

            obj.default_X=parser.Results.default_x;
            if isempty(obj.default_X)
                obj.default_X=zeros(obj.full_len_X,2);
            end
            obj.P=P;
            obj.Wmat=adj;
            obj.adj=adj~=0;
            [obj.I,obj.J]=find(adj);
            od=obj.I~=obj.J;
            obj.I=obj.I(od);
            obj.J=obj.J(od);
            W=full(adj(sub2ind(size(adj),obj.I,obj.J)));
            clamp=1e-2;
            W(W<clamp)=clamp;
            obj.W=W;
%             inds=find(0<obj.Wmat);
%             temp=obj.Wmat(inds)<clamp;
%             inds=inds(temp);
            
%             obj.Wmat(inds)=clamp;
            inds=sub2ind(size(adj),obj.I,obj.J);
            temp=obj.Wmat(inds)<clamp;
            inds=inds(temp);
            obj.Wmat(inds)=clamp;
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
            obj.right_to_left=zeros(length(adj),1)';
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
            
            
            
            obj.fix_gradient_tol= parser.Results.fix_gradient_tol;
            
            %we always fix the vertices which are right side of cut, and
            %cone vertices
            obj.freeVertices=~(obj.right_side_vertices|obj.cones);
            %also fix all vertices requested by user to fix
            obj.freeVertices(parser.Results.fix_inds)=false;
            
            
            
            
           
             if obj.fix_gradient_tol>0
                 %fix also vertices whose gradient is big
                 x=obj.default_X;
                 x=x';
                 x=x(:);
                [~,g]=obj.objective_and_grad(x);
                g=reshape(g,2,round(length(g)/2))';
                n=sqrt(sum(g.^2,2));
                o=obj.karcherPerVertex(x);
                d=n;%./o;
                d(isnan(d))=-inf;
                [~,sortIndex] = sort(d,'descend');  %# Sort the values in
                                                  %#   descending order
                maxIndex = sortIndex(1:1);  %# Get a linear index into A of the 5 largest values
%                 maxIndex=d>obj.fix_gradient_tol;
                obj.freeVertices=false(size(obj.freeVertices));
                obj.freeVertices(maxIndex)=true;
                
                
                
            end
            %             assert(length(obj.Iright)+length(obj.Iinterior)==length(obj.I));
        end
        
        
       
       
        
        
        
        
        
    end
end

