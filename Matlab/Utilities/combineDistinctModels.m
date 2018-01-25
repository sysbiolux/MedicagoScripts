function CombinedModel = combineDistinctModels(model1, model2)
% Simple function to combine two models under the assumption that they have
% nothing in common
combinedFields = intersect(fieldnames(model1),fieldnames(model2));

combinedFields = setdiff(combinedFields,{'S','description','rxnGeneMat','genes','osenseStr','osense'});

CombinedModel = struct();
%Combine all fields
for i = 1:numel(combinedFields)
    fieldName = combinedFields{i};
    CombinedModel.(fieldName) = [model1.(fieldName);model2.(fieldName)];
end
[m1,n1] = size(model1.S);
[m2,n2] = size(model2.S);
%We will ignore the genes for this... 
CombinedModel.S = [model1.S, sparse(m1,n2) ; sparse(m2,n1), model2.S];
CombinedModel.osenseStr = 'max';