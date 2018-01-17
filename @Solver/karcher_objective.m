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