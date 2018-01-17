 function setCutMesh(obj,M_cut)
            if isa(M_cut,'TreeCutter')
                M_cut.V=obj.M_orig.V(:,M_cut.new2old);
                obj.M_cut.V=M_cut.V;
                obj.M_cut.T=M_cut.T;
                if size(obj.M_cut.V,2)~=3
                    obj.M_cut.V=obj.M_cut.V';
                end
                if size(obj.M_cut.T,2)~=3
                    obj.M_cut.T=obj.M_cut.T';
                end
                obj.M_cut.pathPairs=M_cut.pathPairs;
                obj.M_cut.New2Old=M_cut.new2old;
                obj.M_cut.Old2New=M_cut.old2new;
            else
                M_cut.V=obj.M_orig.V(:,M_cut.New2Old);
                obj.M_cut.V=M_cut.V;
                obj.M_cut.T=M_cut.F;
                if size(obj.M_cut.V,2)~=3
                    obj.M_cut.V=obj.M_cut.V';
                end
                if size(obj.M_cut.T,2)~=3
                    obj.M_cut.T=obj.M_cut.T';
                end
                obj.M_cut.pathPairs=M_cut.pathPairs;
                obj.M_cut.New2Old=M_cut.New2Old;
                obj.M_cut.Old2New=M_cut.Old2New;
            end
            
            obj.flat_T=obj.M_cut.T;
            
        end