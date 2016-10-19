function [scans] = ScanDualReactionQP(model,ammoniumname, nitratename ,minv,maxv,steps)
%two reactions (in this instance we assume its nitrate and ammonium) for a range of combined values.
Ammonium  = find(ismember(model.rxns,ammoniumname));
Nitrate = find(ismember(model.rxns,nitratename));
%reac2 = find(ismember(model.rxns,rxnName2));
stepvals1 = linspace(minv,maxv,steps); 
stepvals2 = linspace(maxv,minv,steps); 
for i = 1:steps
    %fix the flux
    model.lb(Ammonium) = stepvals1(i);
    model.ub(Ammonium) = stepvals1(i); 
    model.lb(Nitrate) = stepvals2(i);
    model.ub(Nitrate) = stepvals2(i); 
    lp2 = createQP(model);
    sol = lp2.solve();
    if sol.status == 1
        scans(:,i) = sol.x;
    end
end
end

