classdef PosConstraints<handle
    properties
        A;
        b=[];
        
    end
    
    methods
        function obj=PosConstraints(nvars)
            obj.A=sparse(0,nvars*2);
        end
        function addConstraint(obj,inds,w,rhs)
            assert(length(rhs)==2);
            if size(rhs,1)==1
                rhs=rhs';
            end
            obj.A(end+1,inds*2-1)=w;
            obj.A(end+1,inds*2)=w;
            
            obj.b=[obj.b;rhs];
        end
        function newConstraint(obj)
            obj.A(end+1,1)=0;
            obj.b(end+1)=0;
        end
        function addLineConstraint(obj,ind,n,offset)
            obj.A(end+1,ind*2-1:ind*2)=n;
            
            obj.b=[obj.b;offset];
        end
        
        function addTransConstraints(obj,sinds,tinds,T)
            assert(length(sinds)==length(tinds));
            if sinds(end)==tinds(end)
                sinds=sinds(end:-1:1);
                tinds=tinds(end:-1:1);
            end            
            for ind=2:length(sinds)
                for y=1:2
                    obj.A(end+1,sinds(ind)*2+[-1,0])=T(y,:);
                    obj.A(end,sinds(1)*2+[-1,0])=obj.A(end,sinds(1)*2+[-1,0])-T(y,:);
                    obj.A(end,tinds(ind)*2+y-2)=obj.A(end,tinds(ind)*2+y-2)-1;
                    obj.A(end,tinds(1)*2+y-2)=obj.A(end,tinds(1)*2+y-2)+1;
                    
                    obj.b=[obj.b;0];
                end
            end
        end
        function addMobiusConstraints(obj,sinds,tinds,M,X,side)
            assert(length(sinds)==length(tinds));  
            for ind=2:length(sinds)-1
                if side
                   A=M.a-M.c*(X(tinds(ind),1)+X(tinds(ind),2)*1i);
                   B=-M.d;
               else
                   A=M.a;
                   B=-M.c*(X(sinds(ind),1)+X(sinds(ind),2)*1i)-M.d;
               end
               A=[real(A) -imag(A);
                   imag(A) real(A)];
               B=[real(B) -imag(B);
                  imag(B) real(B)];                
                for y=1:2
                    obj.A(end+1,sinds(ind)*2+[-1,0])=A(y,:);
                    obj.A(end,tinds(ind)*2+[-1,0])=B(y,:);
                end
                obj.b=[obj.b;real(-M.b);imag(-M.b)];
            end
        end
    end
    
end

