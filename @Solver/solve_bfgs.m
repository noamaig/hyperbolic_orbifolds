
        function x=solve_bfgs(obj,X)
            error('deprecated');
            %             X=obj.reduceX(X);
            
            %             opts    = struct( 'factr', 0, 'pgtol', 0, 'm', 10);
            %
            %             opts.x0=X;
            %             [x,~,info]=lbfgsb(@(X)obj.for_bfgs(X),-inf(size(X)),inf(size(X)),opts);
            %              x=obj.expandX(x);
            %             x=obj.correctX(x);
            %             return
           
            options =optimoptions('fmincon');
            options.MaxFunEvals=inf;
            options.Display = 'iter-detailed';
            options.Diagnostics = 'on';
            %             params{SOLVER_NO}.options.DerivativeCheck = 'on';

            options.GradObj = 'on';
            options.MaxIter = 500000;
            options.TolFun=0;1e-10;
            options.TolX=1e-32;
%             if isempty(obj.L)
            
%             else
%                 options.Hessian = 'user-supplied';
%                 options.HessFcn=@(X,lambda)obj.hessFun(X,lambda);
%             end
            %             assert(norm(obj.cones_cons.A*X-obj.cones_cons.b)<1e-10);
            X=obj.reduceX(X);
            if length(X)<10
%                  options.Hessian = 'bfgs';
%                  options.MaxIter = 20;
% %                  options.TolFun=0;
                    options.TolX=1e-16;
            else
                options.Hessian = {'lbfgs',1};
            end
             x = fmincon(@(X)obj.for_bfgs(X),X,[],[],[],[],[],inf(size(X)),[],options);
 
            %             options=[];
            %             options.optTol=1e-10;
            %             options.progTol=1e-12;
            %             options.MaxFunEvals=1e+8;
            %             options.MaxIter=500;
            %             x = minFunc(@(X)obj.for_bfgs(X),X,options);
            x=obj.expandX(x);
            x=obj.correctX(x);
        end