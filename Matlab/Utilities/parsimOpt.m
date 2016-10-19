function sol = parsimOpt(model)
%Do a Parimonous optimization of the model
%Generate the splitted lp
lp = createSplitLP(model);
%Solve it and extract the real solution (i.e. the combination of forward
%and backward reaction
sol = lp.solve();
realsol = sol.x(1:numel(model.rxns)) - sol.x((numel(model.rxns)+1):end);
%Set the bounds for c according to the objective found
model.lb(find(model.c)) = realsol(find(model.c));
model.ub(find(model.c)) = realsol(find(model.c));
%create an lp, that minimized the overall flux.
lp2 = createSplitLP(model);
sol2 = lp2.solve();
%And combine it again.
maxminsol = sol2.x(1:numel(model.rxns)) - sol2.x((numel(model.rxns)+1):end);
sol.x = maxminsol;



