
function x=solve_bfgs_fast(obj,X)

lbfgs_memory = 5;

options =[];
options.MaxIter = 10000;
options.TolFun=0; %1e-10;
options.TolX=1e-16;

X=obj.reduceX(X);
if length(X)<10
    options.TolX=1e-16;
end

X0 = X;

% setup lbfgs solver
fun = @(Y)obj.for_bfgs_fast(Y);
Ab = AffineSpace(zeros(0,length(X0)),zeros(0,1));
x0 = Ab.projectOnto(X0);
% P = Precond_Identity();
P = Precond_Fixed(-Ab.N'*kron(obj.Wmat(obj.freeVertices,obj.freeVertices),eye(2))*Ab.N);
lbfgs_solver = OptimSolverLBFGS_NEW(fun, x0, Ab, P, lbfgs_memory);

% solve
lbfgs_solver.solveTol(options.TolX, options.TolFun, options.MaxIter);
lbfgs_solver.updatePrecond(Precond_Identity());
lbfgs_solver.updateMemory(0);
x = lbfgs_solver.solveTol(options.TolX, options.TolFun, options.MaxIter);

% output
x=obj.expandX(x);
x=obj.correctX(x);
end