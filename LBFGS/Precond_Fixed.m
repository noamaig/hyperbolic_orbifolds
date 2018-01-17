classdef Precond_Fixed < Precond
    
    properties
        P;
    end
    
    methods
        function obj = Precond_Fixed(A)
            t_start = tic;
            fprintf('Init preconditioner... ');
            
            obj.P = SparseLU(A);
            
            obj.elapsed_time = toc(t_start);
            fprintf('Done (%.3g sec)\n', obj.elapsed_time);
        end
        function y = apply(obj, x, z)
            y = obj.P.solve(x);
        end
    end
    
end