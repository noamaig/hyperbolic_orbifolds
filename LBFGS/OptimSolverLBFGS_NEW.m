classdef OptimSolverLBFGS_NEW < handle
    
    properties
        % inputs
        fun;
        x0;
        Ab;
        m;
        P; % precond
        
        % setup
        dim;
        
        % internal state
        x;
        x_f;
        x_fgrad;
        x_prev;
        x_f_prev;
        x_fgrad_prev;
        p;
        t = nan;
        
        % lbfgs memory
        curr_m;
        y_k;
        s_k;
        rao_k;
        alpha_k;
        beta_k;
        gamma_k;
        H_k;
        q;
        r;
        
        % stop condition variables
        stopCntAccept = 2;
        tolXCnt = 0;
        tolFunCnt = 0;
        rel_delta_x;
                
        % log
        f_count = 0; % store function evaluation count
        grad_count = 0; % store gradient evaluation count
        elapsed_time;
        
        % constants
        verbose = 2;
        tol_Ab = 1e-10;
        useLineSearch = true;
        ls_alpha = 0.2;
        ls_beta = 0.5;   
    end
               
    methods
        function obj = OptimSolverLBFGS_NEW(fun,x0,Ab,P,m)
            fprintf('Constructing LBFGS solver\n');
            % copy
            obj.fun = fun;
            obj.x0 = x0;
            obj.Ab = Ab;
            % init precond
            obj.updatePrecond(P);
            % init lbfgs
            obj.initLBFGS(m); 
        end
        
        function report(obj,verbosity,varargin)
            if verbosity<=obj.verbose
                fprintf(varargin{:});
            end
        end
        
        function initLBFGS(obj,m)
            assert(norm(obj.Ab.A*obj.x0-obj.Ab.b,'inf')<=obj.tol_Ab, 'x0 failed to satisfy the linear constraints!');
            obj.x = obj.x0;
            obj.x_prev = inf;
            % setup memory
            obj.updateMemory(m);    
            % evaluate function and gradient
            [obj.x_f, obj.x_fgrad] = obj.funEvaluateValueGrad(obj.x); 
        end
        
        function updateMemory(obj,m)
            obj.m = m;
            obj.curr_m = 0;
            % reset memory
            obj.y_k = zeros(obj.Ab.dim_N, obj.m);
            obj.s_k = zeros(obj.Ab.dim_N, obj.m);
            obj.rao_k = zeros(1, obj.m);
            obj.alpha_k = zeros(1, obj.m);
            obj.beta_k = zeros(1, obj.m);
            % set dummy memory (for H_k initialization)
            obj.y_k(:,1) = 1;
            obj.s_k(:,1) = 1;        
        end
        
        function updatePrecond(obj,P)
            obj.P = P;
        end
        
        function [z_f] = funEvaluateValue(obj,z)
            [z_f] = obj.fun(z);
            obj.f_count = obj.f_count + 1;
        end
        
        function [z_f, z_fgrad] = funEvaluateValueGrad(obj,z)
            [z_f, z_fgrad] = obj.fun(z);
            obj.f_count = obj.f_count + 1;
            obj.grad_count = obj.grad_count + 1;
        end
        
        function iterate(obj)
            % hold previous step
            obj.x_prev = obj.x;
            obj.x_f_prev = obj.x_f;
            obj.x_fgrad_prev = obj.x_fgrad;
           
            %%%%% LBFGS %%%%% (Nocedal&Wright pp 178-179, alg. 7.4+7.5)
            % compute p
            obj.q = obj.Ab.N'*obj.x_fgrad; % map to null-space coordinates
            for ii = 1:obj.curr_m
                obj.alpha_k(ii) = obj.rao_k(ii)*obj.s_k(:,ii)'*obj.q;
                obj.q = obj.q - obj.alpha_k(ii)*obj.y_k(:,ii);
            end
            
            % init r (implicitychoose H_k^0)
            obj.gamma_k = (obj.s_k(:,1)'*obj.y_k(:,1)) / (obj.y_k(:,1)'*obj.y_k(:,1));
            obj.r = obj.gamma_k * obj.P.apply(obj.q, obj.x);
            
            % compute r
            for ii = 1:obj.curr_m
                obj.beta_k(ii) = obj.rao_k(ii)*obj.y_k(:,ii)'*obj.r;
                obj.r = obj.r + obj.s_k(:,ii)*(obj.alpha_k(ii) - obj.beta_k(ii));
            end
            
            % set search direction
            obj.p = -obj.Ab.N*obj.r;  % map back to ambient coordinates

            % line search
            obj.t = 1;
            if obj.useLineSearch
                [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond;
                while linesearch_cond_lhs>linesearch_cond_rhs
                    obj.t = obj.ls_beta*obj.t;
                    [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond;
                end
            end
            
            % update
            obj.x = obj.x + obj.t*obj.p;
            
            % project onto the sphere - hack
            % %             temp = normc(reshape(obj.x,3,[]));
            % %             obj.x = temp(:);
            
            % evaluate function and gradient
            [obj.x_f, obj.x_fgrad] = obj.funEvaluateValueGrad(obj.x);
            
            % update lbfgs memory
            obj.y_k(:,2:end) =  obj.y_k(:,1:end-1);
            obj.s_k(:,2:end) =  obj.s_k(:,1:end-1);
            obj.rao_k(:,2:end) =  obj.rao_k(:,1:end-1);
            obj.curr_m = min(obj.curr_m+1, obj.m);
            obj.y_k(:,1) = obj.Ab.N'*(obj.x_fgrad - obj.x_fgrad_prev);  % map to null-space coordinates
            obj.s_k(:,1) = obj.Ab.N'*(obj.x - obj.x_prev);  % map to null-space coordinates
            obj.rao_k(:,1) = 1/(obj.y_k(:,1)'*obj.s_k(:,1));
                       
            function [linesearch_cond_lhs, linesearch_cond_rhs] = computeLineSearchCond
                linesearch_cond_lhs = obj.funEvaluateValue(obj.x + obj.t*obj.p);
                linesearch_cond_rhs = obj.x_f + obj.ls_alpha*obj.t*obj.x_fgrad'*obj.p;
            end
        end
        
        function [x, x_f, log] = solveTol(obj, TolX, TolFun, max_iter)
            t_start = tic;
            obj.report(1,'Running LBFGS (TolX: %g   TolFun: %g   #max iter: %d   m: %d)\n', TolX, TolFun, max_iter, obj.m);
            printDiagnosticsHeader;
            logInit;
            logState(0);
            for iter = 1:max_iter
                obj.iterate();      
                stop = stopCriteria;
                printDiagnostics;
                logState(iter);
                if stop
                    break;
                end
            end
            if ~stop
                obj.report(1,'Maximal number of iteration exceeded (max_iter)\n', iter);
            end
            % return 
            x = obj.x;
            x_f = obj.x_f;
            logTruncate;
            obj.elapsed_time = toc(t_start);
            fprintf('LBFGS done (%.3g sec)\n', obj.elapsed_time);
            
            function printDiagnosticsHeader
                obj.report(2, '   Iter \t f_cnt \t f(x) \t rel_delta_x \t ||grad(x)||_inf \t gamma_k\n');
            end
                
            function printDiagnostics
                obj.report(2, '   %d \t %d \t %g \t %g \t %g \t %g\n', iter, obj.f_count, obj.x_f, obj.rel_delta_x, norm(obj.x_fgrad,'inf'), obj.gamma_k);
            end
            
            function logInit
                log.iter = nan(1, max_iter+1);
                log.f_count = nan(1, max_iter+1);
                log.f = nan(1, max_iter+1);
                log.t = nan(1, max_iter+1);
                log.exitflag = 'MaxIter';
            end
            function logState(iter)
                log.iter(iter+1) = iter;
                log.f_count(iter+1) = obj.f_count;
                log.f(iter+1) = obj.x_f;
                log.t(iter+1) = obj.t;
            end
            function logTruncate    
                log.iter = log.iter(1:iter+1);
                log.f_count = log.f_count(1:iter+1);
                log.f = log.f(1:iter+1);
                log.t = log.t(1:iter+1);
            end
            
            function stop = stopCriteria
                stop = false;
                if iter>1
                    % No progress at all
                    if norm(obj.x_prev-obj.x,'inf') == 0
                        obj.report(1,'Stopped after %d iterations (NoProgress)\n', iter);
                        log.exitflag = 'NoProgress';
                        stop = true;
                    end          
                    % TolX
                    obj.rel_delta_x = norm(obj.x_prev-obj.x,'fro')/(1+norm(obj.x_prev,'fro'));
                    if obj.rel_delta_x<TolX
                        obj.tolXCnt = obj.tolXCnt + 1;
                    else
                        obj.tolXCnt = 0;
                    end
                    if obj.tolXCnt>=obj.stopCntAccept
                        obj.report(1,'Stopped after %d iterations (TolX)\n', iter);
                        log.exitflag = 'TolX';
                        stop = true;
                    end
                    % TolFun
                    if abs(obj.x_f_prev-obj.x_f)<TolFun*(1+abs(obj.x_f_prev))
                        obj.tolFunCnt = obj.tolFunCnt + 1;
                    else
                        obj.tolFunCnt = 0;
                    end
                    if obj.tolFunCnt>=obj.stopCntAccept
                        obj.report(1,'Stopped after %d iterations (TolFun)\n', iter);
                        log.exitflag = 'TolFun';
                        stop = true;
                    end
                end
            end
        end

    end
    
    
    
end