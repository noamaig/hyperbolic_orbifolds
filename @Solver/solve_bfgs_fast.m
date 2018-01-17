
function x=solve_bfgs_fast(obj,X)

%             X=obj.reduceX(X);

%             opts    = struct( 'factr', 0, 'pgtol', 0, 'm', 10);
%
%             opts.x0=X;
%             [x,~,info]=lbfgsb(@(X)obj.for_bfgs(X),-inf(size(X)),inf(size(X)),opts);
%              x=obj.expandX(x);
%             x=obj.correctX(x);
%             return

lbfgs_memory = 5;

options =[];%optimoptions('fmincon');
options.MaxFunEvals=inf;
options.Display = 'iter-detailed';
options.Diagnostics = 'on';
%             params{SOLVER_NO}.options.DerivativeCheck = 'on';

options.GradObj = 'on';
options.MaxIter = 500000;
options.TolFun=0;1e-10;
options.TolX=1e-16;
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


X0 = X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup lbfgs solver
t_start = tic;
fun = @(Y)obj.for_bfgs_fast(Y);
Ab = AffineSpace(zeros(0,length(X0)),zeros(0,1));
x0 = Ab.projectOnto(X0);
% P = Precond_Identity();
P = Precond_Fixed(-Ab.N'*kron(obj.Wmat(obj.freeVertices,obj.freeVertices),eye(2))*Ab.N);
lbfgs_solver = OptimSolverLBFGS_NEW(fun, x0, Ab, P, lbfgs_memory);
log.t_init_lbfgs = toc(t_start);


% solve
t_start = tic;
[x,~,log_lbfgs] = lbfgs_solver.solveTol(options.TolX, options.TolFun, options.MaxIter);
log.t_lbfgs = toc(t_start);
log.exitflag = log_lbfgs.exitflag;
log.n_iter = log_lbfgs.iter(end);
log.f_count = log_lbfgs.f_count(end);

% plot some stuff
if false
    figure;
    plot(log_lbfgs.iter,log_lbfgs.f);
    ylabel('f');
    
    figure;
    plot(log_lbfgs.iter,log_lbfgs.f_count);
    ylabel('f_count');
    
    figure;
    plot(log_lbfgs.iter,log_lbfgs.t);
    ylabel('t');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% compare to fmincon %%
% x_old = fmincon(@(X)obj.for_bfgs_fast(X),X0,[],[],[],[],[],inf(size(X)),[],options);
% norm(x_old-x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


%             options=[];
%             options.optTol=1e-10;
%             options.progTol=1e-12;
%             options.MaxFunEvals=1e+8;
%             options.MaxIter=500;
%             x = minFunc(@(X)obj.for_bfgs(X),X,options);
x=obj.expandX(x);
x=obj.correctX(x);
end