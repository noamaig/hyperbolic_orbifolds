classdef Precond_Identity < Precond
    
    properties
    end
    
    methods
        function y = apply(obj, x, z)
            y = x;
        end
    end
    
end