classdef AffineSpace < handle
    
    properties
        A;
        b;
        N;
        dim_N;
        x_p;
        
        tol = 1e-8;
        
        elapsed_time;
    end
    
    
    methods
        function obj = AffineSpace(A,b)
            t_start = tic;
            fprintf('Computing affine space... ');
            
            obj.A = A;
            obj.b = b;
            
            if ~isempty(A)
                % obj.N = sparse(null(full(obj.A)));
                % warning('Verify this null-space computation / or replace...');
                %             assert(svds(A,1,'sm')>obj.tol, 'A must be full rank!');
                [Q,~] = qr(A');
                obj.N = Q(:,size(A,1)+1:end);
                
                obj.dim_N = size(obj.N,2);
                obj.x_p = A\b;
            else
                obj.N = speye(size(A,2));
                obj.dim_N = size(A,2);
                obj.x_p = 0;
            end
            
            obj.elapsed_time = toc(t_start);
            fprintf('Done (%.3g sec)\n', obj.elapsed_time);
        end
        
        function y = projectOnto(obj,x)
            y = obj.x_p + obj.N*(obj.N'*(x-obj.x_p));
        end
    end
    
    
end