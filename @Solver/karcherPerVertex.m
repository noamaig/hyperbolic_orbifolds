
        function o=karcherPerVertex(obj,x)
            persistent I1 J1 I2 J2 inds E;
        
        if size(x,2)==1
                x=[x(1:2:end) x(2:2:end)];
            end
            if isempty(E)
            %all edges
            [I,J]=find(obj.adj);
            %remove the diagonal
            inds=I~=J;
            I1=I(inds);
            J1=J(inds);
            
            inds=obj.right_to_left(I1)~=0;
            I2=I1;
            J2=J1;
            I2(inds)=obj.right_to_left(I1(inds));
            inds=obj.right_to_left(J1)~=0;
            J2(inds)=obj.right_to_left(J1(inds));
            
            E=[I2 J2];
            [E,inds]=unique(E,'rows');
            end
            %take the weights
            w=obj.Wmat(sub2ind(size(obj.Wmat),I1,J1));
            %compute the distance for each edge
            a=x(I1,1);
            b=x(I1,2);
            c=x(J1,1);
            d=x(J1,2);
            [~,~,~,~,O]=obj.karcher_grad(a,b,c,d,w);
            %the weighted squared distance
            O=w.*O.^2;
            %convert to edges of original uncut mesh
            
            
            O=O(inds);
            I3=I2(inds);
            o=accumarray(I3,O);
            
            o(obj.right_to_left~=0)=o(obj.right_to_left(obj.right_to_left~=0));
        end