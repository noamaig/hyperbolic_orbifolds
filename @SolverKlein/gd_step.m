function [newY,y_f]=gd_step(obj,y)
% line search
beta=0.1;
alpha=1e-2;
[y_f,y_fgrad]=obj.for_bfgs(y);
p=-y_fgrad;
t = 1;
[linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond;
while linesearch_cond_lhs>linesearch_cond_rhs 

    t = beta*t;
    [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond;
end

newY=y+t*p;
    function [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond
        linesearch_cond_lhs =obj.for_bfgs(y + t*p);
        linesearch_cond_rhs = y_f + alpha*t*y_fgrad'*p;
    end
end