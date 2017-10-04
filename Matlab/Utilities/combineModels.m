function CombinedModel = combineModels(model1, model2)
% Simple function to combine two models under the assumption, 
% that they share the set of genes And they have no metabolites in common
combinedFields = intersect(fieldnames(model1),fieldnames(model2));

combinedFields = setdiff(combinedFields,{'S','rxnGeneMat','genes'});

CombinedModel = struct();
%Combine all fields
for i = 1:numel(combinedFields)
    fieldName = combinedFields{i};
    CombinedModel.(fieldName) = [model1.(fieldName);model2.(fieldName)];
end
[m1,n1] = size(model1.S);
[m2,n2] = size(model2.S);
if isfield(model1,'rxnGeneMat') && isfield(model2,'rxnGeneMat')
    CombinedModel.rxnGeneMat = [model1.rxnGeneMat;model2.rxnGeneMat];
end
CombinedModel.S = [model1.S, sparse(m1,n2) ; sparse(m2,n1), model2.S];
CombinedModel.genes = model1.genes;