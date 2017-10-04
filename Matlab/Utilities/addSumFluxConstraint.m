function model = addSumFluxConstraint(model,reactions,coefs,lb,ub,name)
%add a sumfluxconstraint. This will introduce the blb and bub fields to the
%model. i.,e. lb <= sum(flux_through reactions) <=  ub
if ~isfield(model, 'blb')
    model.blb = model.b;
    model.bub = model.b;
end

reacpos = find(ismember(model.rxns,reactions));
model = addMetabolite(model,name,name,'','','','','',0,0);
%This model is no longer useable by optimizecbmodel...
model.bub(end) = ub;
model.blb(end) = lb;
%Add the metabolite to all reactions.
model.S(end,reacpos) = coefs;

