classdef SolverKlein <handle
    properties
        inds;
        W;
       
        I;
        J;
       
        
        adj;
      
        
        L;
        fix_gradient_tol;
        freeVertices;
        full_len_X;
        default_X;
        Wmat;
        precondition=false;
        radius=1;
        A;
        b;
    end
    methods
        
        
       
        
        
       
        
        
        
      
        function obj=SolverKlein(adj,boundary_paths,boundary_normals,boundary_shifts)
            parser = inputParser;
            parser.addOptional('fix_gradient_tol',0,@isnumeric);
            parser.addOptional('fix_inds',[],@(X)all(isnumeric(X)));
            parser.addOptional('default_x',zeros(length(adj),2));
            parser.parse([]);
                         obj.full_len_X=length(adj);

            obj.default_X=parser.Results.default_x;
            if isempty(obj.default_X)
                obj.default_X=zeros(obj.full_len_X,2);
            end
%             obj.boundary_paths=boundary_paths;
%             obj.boundary_normals=boundary_normals;
%             obj.boundary_shifts=boundary_shifts;
            obj.Wmat=adj;
            obj.adj=adj~=0;
            [obj.I,obj.J]=find(adj);
            od=obj.I~=obj.J;
            obj.I=obj.I(od);
            obj.J=obj.J(od);
            W=full(adj(sub2ind(size(adj),obj.I,obj.J)));
            clamp=1e-1;
            W(W<clamp)=clamp;
            obj.W=W;
%             warning('add half weights on boundary!');
            inds=find(0<obj.Wmat);
            temp=obj.Wmat(inds)<clamp;
            inds=inds(temp);
            
            obj.Wmat(inds)=clamp;
            
            %we always fix the vertices which are right side of cut, and
            %cone vertices
            obj.freeVertices=true(length(adj),1);
            
            %want to enforce V(b,:)*n=a
            J=[];
            V=[];
            I=[];
            b=[];
            Neq=0;
            for i=1:length(boundary_paths)
                p=boundary_paths{i};
                J1=p*2-1;
                V1=boundary_normals(i,1);
                J2=p*2;
                V2=boundary_normals(i,2);
                J=[J;J1;J2];
                V=[V;repmat(V1,length(J1),1);repmat(V2,length(J2),1)];
                I=[I;Neq+repmat(1:length(p),1,2)'];
                Neq=max(I);
                b=[b;repmat(boundary_shifts(i),length(p),1)];
                
            end
            obj.A=sparse(I,J,V,max(I),length(adj)*2);
            obj.b=b;
            assert(length(obj.b)==size(obj.A,1));
        end
        
        
       
       
        
        
        
        
        
    end
end

