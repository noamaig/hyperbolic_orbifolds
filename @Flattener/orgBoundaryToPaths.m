function orgBoundaryToPaths(obj,N)
            pathPairs={};
            inds=obj.orgBoundary;
            d=sqrt(sum(obj.M_orig.V(:,inds(1:end-1))'-obj.M_orig.V(:,inds(2:end))',2).^2);
            d=[0;cumsum(d)/sum(d)];
            
            inds=[1];
            for i=1:N-1
                ind=find(d>i/N,1);
                inds(end+1)=ind;
            end
            %             inds(end+1)=length(obj.orgBoundary);
            
            for i=1:length(inds)-1
                p=obj.orgBoundary(inds(i):inds(i+1));
                pathPairs{i}=[p p];
            end
            p=obj.orgBoundary([inds(i+1):length(obj.orgBoundary) 1]);
            pathPairs{end+1}=[p p];
            obj.uncut_cone_inds=obj.orgBoundary(inds);
            M_cut=[];
            M_cut.V=obj.M_orig.V';
            M_cut.T=obj.M_orig.F';
            M_cut.pathPairs=pathPairs;
            M_cut.Old2New=num2cell(1:length(obj.M_orig.V));
            M_cut.New2Old=1:length(obj.M_orig.V);
            obj.M_cut=M_cut;
        end
       