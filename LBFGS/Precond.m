classdef (Abstract) Precond < handle
    
    properties
        elapsed_time;
    end
    
    methods (Abstract)
        y = apply(obj, x, z) % return y = P_z (x), where P is a preconditioner acting on x, which depends on z
    end
    
end