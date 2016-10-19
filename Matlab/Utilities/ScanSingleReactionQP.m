function [scans] = ScanSingleReactionQP(model,rxnName,minv,maxv,steps)
%Scan Reaction using QP
reac  = find(ismember(model.rxns,rxnName));
%reac2 = find(ismember(model.rxns,rxnName2));
stepvals = linspace(minv,maxv,steps); 
for i = 1:steps
    %fix the flux
    model.lb(reac) = stepvals(i);
    model.ub(reac) = stepvals(i); 
    lp2 = createQP(model);
    sol = lp2.solve();
    if sol.status == 1
        scans(:,i) = sol.x;
    end
end
end

