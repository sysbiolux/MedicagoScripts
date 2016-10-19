function splitLP = createSplitLP(model)
%This function assumes that there are two additional fields in the model
%structure:
%bub and blb for upper and lower constraint bounds.

splitLP = Cplex();
A = [model.S, -model.S];
[m,n] = size(A);
%lbs: 
%zero or -model.lb for all forward reactions.
%zero or -model.ub

lbs = [max(0,model.lb);max(0,-model.ub)];

%ubs
ubs = [max(0,model.ub) ; max(0,-model.lb)];
if isfield(model,'bub')
    constraintpos = find(model.bub);
    A(constraintpos,:) = abs(A(constraintpos,:));
    rhs = model.bub;
    lhs = model.blb;
else
    rhs = zeros(size(A,1),1);
    lhs = zeros(size(A,1),1);
end
obj = [model.c ; -model.c];
splitLP.Model.sense = 'maximize';
%rownames = model.mets;
%colnames = [model.rxns ; strcat(model.rxns,'_r')];
splitLP.Model.A = A;
splitLP.Model.obj = obj;
splitLP.Model.lb = lbs;
splitLP.Model.ub = ubs;
splitLP.Model.lhs = lhs;
splitLP.Model.rhs = rhs;
%splitLP.Model.rowname = rownames;
%splitLP.Model.colname = colnames;
splitLP.DisplayFunc = [];
end
