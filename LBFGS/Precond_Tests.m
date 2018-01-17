classdef Precond_Tests < Precond
    
    properties
        N;
        T;
        V;
        meshDiff;
        origAreas;
        currAreas;
        L;
    end
    
    methods
        function obj = Precond_Tests(N,T,V)
            t_start = tic;
            fprintf('Init preconditioner... ');
            obj.N = N;
            obj.T = T;
            obj.V = V;
            
            [obj.meshDiff, obj.origAreas] = getFlatteningDiffCoefMatrix(obj.V,obj.T);
            obj.meshDiff = obj.meshDiff(:,1:end/2);
            
            
            obj.elapsed_time = toc(t_start);
            fprintf('Done (%.3g sec)\n', obj.elapsed_time);
        end
        function y = apply(obj, x, z)
            [~, obj.currAreas] = getFlatteningDiffCoefMatrix(reshape(z,3,[])',obj.T);
            
            
            
            
% % %             [TT,areas] = getFlatteningDiffCoefMatrix(V,T);
% % %             TT = TT(:,1:end/2);
            a = obj.origAreas./(obj.currAreas/mean(obj.currAreas));
            wTT = kron(sqrt(a), ones(2^2,1)).*obj.meshDiff;
            obj.L = -(wTT'*wTT);
            
            figure(300); cla;
            plot(a);
            drawnow;
            
            
            %obj.L = cotmatrix(obj.V,obj.T);
            %obj.L = -obj.meshDiff'*obj.meshDiff;
            y = (-obj.N'*kron(obj.L ,eye(3))*obj.N)\x;
        end
    end
    
end