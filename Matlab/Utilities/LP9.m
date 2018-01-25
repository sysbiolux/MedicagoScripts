function V = LP9(K, P, model, epsilon)
% CPLEX implementation of LP-9 for input sets K, P (see FASTCORE paper)
%
% USAGE:
%
%    V = LP9(K, P, model, epsilon)
%
% .. Authors: -  Nikos Vlassis, Maria Pires Pacheco, Thomas Sauter, 2013
%             LCSB / LSRU, University of Luxembourg

    scalingfactor = 1e5;
    oldTol = getCobraSolverParams('LP','feasTol');
    changeCobraSolverParams('LP','feasTol',oldTol*scalingfactor);
    finishUp = onCleanup(@() changeCobraSolverParams('LP','feasTol',oldTol));
    
    V = [];
    if isempty(P) || isempty(K)
        return;
    end

    np = numel(P);
    nk = numel(K);
    [m,n] = size(model.S);

    % objective
    f = [zeros(n,1); ones(np,1)];

    % equalities
    Aeq = [model.S, sparse(m,np)];
    beq = zeros(m,1);

    % inequalities
    Ip = sparse(np,n); Ip(sub2ind(size(Ip),(1:np)',P(:))) = 1;
    Ik = sparse(nk,n); Ik(sub2ind(size(Ik),(1:nk)',K(:))) = 1;
    Aineq = sparse([[Ip, -speye(np)]; ...
                    [-Ip, -speye(np)]; ...
                    [-Ik, sparse(nk,np)]]);
    bineq = [zeros(2*np,1); -ones(nk,1)*epsilon*scalingfactor];

    % bounds
    lb = [model.lb; zeros(np,1)] * scalingfactor;
    ub = [model.ub; max(abs(model.ub(P)),abs(model.lb(P)))] * scalingfactor;

    % Set up LP problem
    LPproblem.A=[Aeq;Aineq];
    LPproblem.b=[beq;bineq];
    LPproblem.lb=lb;
    LPproblem.ub=ub;
    LPproblem.c=f;
    LPproblem.osense=1;%minimise
    LPproblem.csense(1:size(LPproblem.A,1))='E';
    LPproblem.csense(size(Aeq,1)+1:size(LPproblem.A,1))='L';

    solution = solveCobraLP(LPproblem);

    if solution.stat~=1
        fprintf('\n%s%s\n',num2str(solution.stat),' = sol.stat')
        fprintf('%s%s\n',num2str(solution.origStat),' = sol.origStat')
        warning('LP solution may not be optimal')
    end

    x=solution.full;

    if ~isempty(x)
        V = x(1:n);
    else
        V=ones(n,1)*NaN;
    end
    