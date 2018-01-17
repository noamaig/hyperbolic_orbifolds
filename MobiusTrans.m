classdef MobiusTrans <handle
    %% Object representing a Mobius transformation
    
    properties
        a;
        b;
        c;
        d;
    end
    
    methods
        function obj=MobiusTrans(a,b,c,d)
            %constructor - construct the mobius trans f(z)=(a*z+b)/(c*z+d)
            obj.a=a./d;
            obj.b=b./d;
            obj.c=c./d;
            obj.d=1;
        end
        function M2=inverse(obj)
            %construct the inverse transformation
            M2=MobiusTrans(obj.d,-obj.b,-obj.c,obj.a);
        end
%         function A=subsref(obj,S)
%             assert(length(S)==1);
%             assert(strcmp(S.type,'()'));
%             inds=S.subs{1};
%             
%             A=MobiusTrans(obj.a(inds),obj.b(inds),obj.c(inds),obj.d(inds));
%         end
        function M2=compose(other,obj)
            %compose this mobius transformation on the other transfroatmion
            a=other.a.*obj.a+other.b.*obj.c;
            b=other.a.*obj.b+other.b.*obj.d;
            c=other.c.*obj.a+other.d.*obj.c;
            d=other.c.*obj.b+other.d.*obj.d;
            M2=MobiusTrans(a,b,c,d);
        end
        function q=map(obj,p)
            %map points (represented as array of complex numbers)
            was2d=false;
            if size(p,2)==2
                was2d=true;
                p=p(:,1)+1i*p(:,2);
            end
            q=(obj.a.*p+obj.b)./(obj.c.*p+obj.d);
            if was2d
                q=[real(q) imag(q)];
            end
        end
        function d=diff(obj,other)
            %difference (in coefficients) between two mobius trans
            d= max(abs([obj.a-other.a obj.b-other.b obj.c-other.c obj.d-other.d]));
        end
    
    end
    
end

