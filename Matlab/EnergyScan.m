function results = EnergyScan()
% This function replicates the energy usage results as indicated in Supplemental Material 2 
scriptPath = fileparts(which(mfilename));
origDir = cd(scriptPath);
addpath([scriptPath filesep 'Utilities']);

medicago = importMedicago();

carbon_importers = {'TEC_CARBON-DIOXIDE','TEC_GLC', 'TEC_SUCROSE','TEH_Light','TEH_Starch'};
% Close all energy uptake
closeEnergyModel = changeRxnBounds(medicago,carbon_importers, 0, 'u');

% Open Starch import
openStarch = changeRxnBounds(closeEnergyModel,'TEH_Starch', 99999, 'u');

% optimise starch consumption
openStarch = changeObjective(openStarch, 'TEH_Starch',1);

% Open DEHIG AND ATPase flux (free redox equivalents and energy).
% But do not allow them to be used to waste energy (that should be possible
% in the model anyways)
atpasepos = find(ismember(openStarch.rxns,'ATPase'));
openStarch.lb(atpasepos) = -99999;
openStarch.ub(atpasepos) = 99999;
dehogpos = find(ismember(openStarch.rxns,'DEHOG'));
openStarch.lb(dehogpos) = -99999;
openStarch.ub(dehogpos) = 99999;



% get exporters:
pattern = '^(THE_|TCE_|TGE_)';
exchangers = findReactionsWithRegexp(openStarch,pattern);

ammoniumModel = changeRxnBounds(changeRxnBounds(openStarch,'TEC_NITRATE',0,'b'),'TEC_AMMONIUM',99999,'u');
nitrateModel = changeRxnBounds(changeRxnBounds(openStarch,'TEC_AMMONIUM',0,'b'),'TEC_NITRATE',99999,'u');

result = {};
for i=1:numel(exchangers)
    met = findMetsFromRxns(medicago,exchangers{i});
    met = met{1};
    metName = medicago.metNames(find(ismember(medicago.mets,met)));
    currentModel = changeRxnBounds(ammoniumModel,exchangers{i},1,'b');
    res = calcATPAndNADHMin(currentModel);
    % Clean up results for numeric issues
    res.v(abs(res.v)<1e-8) = 0;
    currentModel2 = changeRxnBounds(nitrateModel,exchangers{i},1,'b');
    res2 = calcATPAndNADHMin(currentModel2);
    res2.v(abs(res2.v)<1e-8) = 0;
    
    element = {met,metName,-res.v(dehogpos),-res.v(atpasepos),-res2.v(dehogpos),-res2.v(atpasepos)};
    result{end+1} = element;
end
  
results = cell2table(vertcat(result{:}), "VariableNames", ["Metabolite","Common Name","NADH_Ammonium","ATP_ammonium","NADH_Nitrate","ATP_Nitrate"]);

rmpath([scriptPath filesep 'Utilities']);
cd(origDir)
end 



