function computeV2A(obj)
            fprintf('==== Computing V2A matrix ===\n');
            tid=tic;
            [obj.V2A,obj.areas] = getFlatteningDiffCoefMatrix(obj.M_cut.V,obj.M_cut.T); % calculate map between 2d vertices to differentials
            fprintf('done computing, %f seconds\n',toc(tid));
        end