function sol = parsimOptQP(model);
%Get a solution for a problem (defined in a COBRA model structure) with 
%a minimal quadratic flux sum.

lp = createSplitLP(model);
sol = lp.solve();
realsol = sol.x(1:numel(model.rxns)) - sol.x((numel(model.rxns)+1):end);
sol1 = sol;
model.lb(find(model.c)) = realsol(find(model.c));
model.ub(find(model.c)) = realsol(find(model.c));
lp2 = createQP(model);
sol = lp2.solve();
sol.objval = sol1.objval; 



