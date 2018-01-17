classdef Precond_TargetLaplacian < Precond
    
    properties
        N;
        T;
        L;
    end
    
    methods
        function obj = Precond_TargetLaplacian(N,T)
            t_start = tic;
            fprintf('Init preconditioner... ');
            obj.N = N;
            obj.T = T;
            
            obj.elapsed_time = toc(t_start);
            fprintf('Done (%.3g sec)\n', obj.elapsed_time);
        end
        function y = apply(obj, x, z)
            obj.L = cotmatrix(reshape(z,3,[])',obj.T);
            y = (-obj.N'*kron(obj.L ,eye(3))*obj.N)\x;
        end
    end
    
end