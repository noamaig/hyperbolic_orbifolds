classdef UncutSurfMap <handle
    % a mapping of the vertices of the uncut meshes
    
    properties
        barCoords;
        flatteners;
    end
    
    methods
        function obj=UncutSurfMap(flatteners)
            
            % BC is for vertices on cut mesh. We care only about original mesh. So we can
            % take the BC of only one of the duplucated vertices for each of the old
            % vertices.
            obj.flatteners=flatteners;
            obj.barCoords=cell(length(flatteners));
        end
        function compute(obj,I,J)
           %compute mapping of I --> J
            flattener1=obj.flatteners{I};
            flattener2=obj.flatteners{J};
            M1_uncut=flattener1.M_orig;
            M2_uncut=flattener2.M_orig;
            
            
            
            
            M1_cut=flattener1.M_cut;
            M2_cut=flattener2.M_cut;
            O2N1=M1_cut.Old2New;
            N2O2=M2_cut.New2Old;
            O2N1=cellfun(@(x)x(end),O2N1); % get the first entry of each cell
            [~,~,BC1to2]=flattener2.lift(flattener1);%for each point on 1, barcoords on 2
            temp1=BC1to2(O2N1,:);
            
            
            M=sparse(size(M1_uncut.V,2),size(M2_uncut.V,2));
            
            for i=1:size(temp1,1)
                ti=temp1(i,:);
                ind=find(ti);
                M(i,N2O2(ind))=ti(ind);
            end
            obj.barCoords{J,I}=M;
            
        end
        
            
        
    end
end

