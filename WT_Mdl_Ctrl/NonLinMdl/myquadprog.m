function [X,fval] = myquadprog(H,f,A,B,varargin)
%QUADPROG Quadratic programming.
%   X = QUADPROG(H,f,A,b) attempts to solve the quadratic programming
%   problem:
%
%            min 0.5*x'*H*x + f'*x   subject to:  A*x <= b
%             x
%
%   X = QUADPROG(H,f,A,b,Aeq,beq) solves the problem above while
%   additionally satisfying the equality constraints Aeq*x = beq. (Set A=[]
%   and B=[] if no inequalities exist.)
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that the solution is in the
%   range LB <= X <= UB. Use empty matrices for LB and UB if no bounds
%   exist. Set LB(i) = -Inf if X(i) is unbounded below; set UB(i) = Inf if
%   X(i) is unbounded above.
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0) sets the starting point to X0.
%
%   X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS) minimizes with the
%   default optimization parameters replaced by values in OPTIONS, an
%   argument created with the OPTIMOPTIONS function. See OPTIMOPTIONS for
%   details.
%
%   X = QUADPROG(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with matrix 'H' in PROBLEM.H, the vector 'f' in PROBLEM.f,
%   the linear inequality constraints in PROBLEM.Aineq and PROBLEM.bineq,
%   the linear equality constraints in PROBLEM.Aeq and PROBLEM.beq, the
%   lower bounds in PROBLEM.lb, the upper bounds in PROBLEM.ub, the start
%   point in PROBLEM.x0, the options structure in PROBLEM.options, and
%   solver name 'quadprog' in PROBLEM.solver. Use this syntax to solve at
%   the command line a problem exported from OPTIMTOOL.
%
%   [X,FVAL] = QUADPROG(H,f,A,b) returns the value of the objective
%   function at X: FVAL = 0.5*X'*H*X + f'*X.
%
%   [X,FVAL,EXITFLAG] = QUADPROG(H,f,A,b) returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%
%   All algorithms:
%     1  First order optimality conditions satisfied.
%     0  Maximum number of iterations exceeded.
%    -2  No feasible point found.
%    -3  Problem is unbounded.
%   Interior-point-convex only:
%     2  Solver stalled at feasible point.
%    -6  Non-convex problem detected.
%    -8  Unable to compute step direction; no further progress can be made.
%   Trust-region-reflective only:
%     3  Change in objective function too small.
%    -4  Current search direction is not a descent direction; no further
%        progress can be made.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = QUADPROG(H,f,A,b) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations,
%   maximum of constraint violations in OUTPUT.constrviolation, the
%   type of algorithm used in OUTPUT.algorithm, the number of conjugate
%   gradient iterations (if used) in OUTPUT.cgiterations, a measure of
%   first order optimality (large-scale algorithm only) in
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = QUADPROG(H,f,A,b) returns the set of
%   Lagrangian multipliers LAMBDA, at the solution: LAMBDA.ineqlin for the
%   linear inequalities A, LAMBDA.eqlin for the linear equalities Aeq,
%   LAMBDA.lower for LB, and LAMBDA.upper for UB.
%
%   See also LINPROG, LSQLIN.

%   Copyright 1990-2019 The MathWorks, Inc.

options = struct( ...
    'Algorithm','interior-point-convex', ...
    'Diagnostics','off', ...
    'Display','off', ...
    'HessMult',[], ...
    'MaxIter',200, ... []
    'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
    'PrecondBandWidth',0, ...
    'ProblemdefOptions', struct, ...
    'TolCon',1e-8, ...
    'TolFun',1e-8,...%[], ...
    'TolFunValue', [], ...
    'TolPCG',0.1, ...
    'TolX',1e-12, ...100*eps, ...
    'TypicalX','ones(numberOfVariables,1)', ...
    'LinearSolver', 'auto' ...
    );

X0 = [];

computeFirstOrderOpt = false;
detailedExitMsg = false; %contains(display,'detailed');
verbosity = 0;


numberOfVariables = length(f);
Aeq = zeros(0,numberOfVariables);
Beq = zeros(0,1);

% [X0,lb,ub] = checkbounds(X0,lb,ub,numberOfVariables);
lb = -inf(numberOfVariables,1);
ub = inf(numberOfVariables,1);
defaultopt.ConvexCheck = 1;

% If the output structure is requested, we must reconstruct the
% Lagrange multipliers in the postsolve. Therefore, set computeLambda
% to true if the output structure is requested.
flags.computeLambda = computeFirstOrderOpt;
flags.detailedExitMsg = detailedExitMsg;
flags.verbosity = verbosity;
flags.caller = 'quadprog';

% Full QP
%     if (autoSelect && ~issparse(H)) || strcmp(linearSolver, 'dense')
[X, fval] = ...
    ipqpdense(full(H), f, A, B, Aeq, Beq, lb, ub, X0, flags, ...
    options, defaultopt);

