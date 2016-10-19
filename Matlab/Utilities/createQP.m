function splitLP = createQP(model)
%This function assumes that there are two additional fields in the model
%structure:
%bub and blb for upper and lower constraint bounds.

splitLP = Cplex();
A = [model.S];
[m,n] = size(A);
%lbs: 
%zero or -model.lb for all forward reactions.
%zero or -model.ub

lbs = model.lb;
ubs = model.ub;
if isfield(model,'bub')    
    rhs = model.bub;
    lhs = model.blb;
else
    rhs = zeros(size(A,1),1);
    lhs = zeros(size(A,1),1);
end
%The Objective is set to sum(c_i * (x_i)^2) with c_i >= 0
obj = zeros(n,1);
Q = sparse(diag(model.c));

splitLP.Model.sense = 'minimize';
splitLP.Model.A = A;
splitLP.Model.obj = obj;
splitLP.Model.lb = lbs;
splitLP.Model.ub = ubs;
splitLP.Model.lhs = lhs;
splitLP.Model.rhs = rhs;
splitLP.Model.Q = Q;
splitLP.Param.qpmethod.Cur = 2;
%Turn of the display function
splitLP.DisplayFunc = [];
end
