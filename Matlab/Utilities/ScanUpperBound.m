function [scans] = ScanUpperBound(model,rxnName,minv,maxv,steps)
%Scan a reaction from a minimal upper bound value to a maximal value using linear
%programming,

reac  = find(ismember(model.rxns,rxnName));
stepvals = linspace(minv,maxv,steps); 
for i = 1:steps
    
    model.ub(reac) = stepvals(i);    
    sol = parsimOpt(model);    
    if sol.status == 1
        scans(:,i) = sol.x;
    end
end
end