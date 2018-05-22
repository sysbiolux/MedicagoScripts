function [SymbioticModel,AlaAdjustedTissueModel] = BuildSymbioticModel(TissueModel)
%% uild the Symbiotic Model from the TissueModel (or from the XML Files)
%REad the SMeliloti xml file
%readCbModel needs to initialize tmpCharge!!
Smel = readCbModel(['Data' filesep 'Smel.xml']);
%The following differs for newer Cobra Versions, as there are differences
%in SBML IO (newer versions are more stringent. 

scriptPath = fileparts(which(mfilename));
origDir = cd(scriptPath);
addpath([scriptPath filesep 'Utilities']);


%Replace the Reaction Names and Metabolite Names by something "readable"
Smel.rxns = strrep(Smel.rxnNames,' S','_S');
Smel.mets = regexprep(Smel.metNames,'(_S$)','[S]');
Smel.mets = regexprep(Smel.mets,'(_I$)','[I]');
%There are some compounds in the intermembrane Space
%Smel.mets{975} = 'PROTON[I]';
OrigSmel = Smel;
%% Add Exchange Reactions
%Add exchangers To test the Model
Uptake = {'NITROGEN-MOLECULE','OXYGEN-MOLECULE','SULFATE','Pi','CARBON-DIOXIDE','MG+2','FE+2','CPD-3','CO+2','WATER','PROTON','OXALACETIC_ACID','MAL'};
Export = {'WATER','OXYGEN-MOLECULE','CARBON-DIOXIDE','PROTON','HYDROGEN-MOLECULE'};

for i=1:numel(Export)
    Smel = addReaction(Smel,['TSE_' Export{i}],{[ Export{i} '[S]'] }, [-1], 0,0,1000);
end
for i=1:numel(Uptake)
    Smel = addReaction(Smel,['TES_' Uptake{i}],{[ Uptake{i} '[S]'] }, [ 1], 0,0,1000);
end
%save('Smel.mat','Smel','OrigSmel');

%% Link the Model to the Mtrunctula Root

if nargin == 0
    disp('Building Medicago Model');    
    CombinedModel = BuildTissueModel();
    
else
    CombinedModel = TissueModel;
end
Orig_Combined_Model = CombinedModel;
disp('Generating Exchangers between Symbiont and Plant')
%Define the exchanges between Plant and Symbiont
RootToRhizo = {'VAL','LEU','ILE','WATER','CARBON-DIOXIDE','OXYGEN-MOLECULE'};

RhizoToRoot = {'AMMONIA','CARBON-DIOXIDE','WATER','L-ALPHA-ALANINE'};

SymbioModel = combineDistinctModels(CombinedModel,OrigSmel);



for i=1:numel(RhizoToRoot)
    SymbioModel = addReaction(SymbioModel,['TSyR_' RhizoToRoot{i}],{[ RhizoToRoot{i} '[S]'], ['Root_' RhizoToRoot{i} '[C]']}, [-1 1], 0,0,1000);
end
for i=1:numel(RootToRhizo)
    SymbioModel = addReaction(SymbioModel,['TRSy_' RootToRhizo{i}],{['Root_' RootToRhizo{i} '[C]'], [ RootToRhizo{i} '[S]'] }, [-1 1], 0,0,1000);
end

OxygenPos = find(ismember(SymbioModel.rxns,'TRSy_OXYGEN-MOLECULE')); 
%Add additional Transporters (which need to be charge balanced:
%Glutamate from Root to Rhizo
SymbioModel = addReaction(SymbioModel,'TRSy_GLT',{'Root_GLT[C]','Root_PROTON[C]','GLT[S]','PROTON[S]'},[-1 -1 1 1],1,0,1000);
%Aspartate from Symbiont to Root
SymbioModel = addReaction(SymbioModel,'TSyR_L-ASPARTATE',{'Root_L-ASPARTATE[C]','Root_PROTON[C]','L-ASPARTATE[S]','PROTON[S]'},[1 1 -1 -1],1,-1000,1000);
%Ammonium from intermembrane space:
%SymbioModel = addReaction(SymbioModel,'TIR_AMMONIA',{'Root_AMMONIA[C]','Root_PROTON[C]','AMMONIA[I]','PROTON[I]'},[1 1 -1 -1],0,0,1000);
%Direct Ammonium transport
SymbioModel = addReaction(SymbioModel,'TSyR_AMMONIUM',{'Root_AMMONIA[C]','Root_PROTON[C]','AMMONIUM[S]','PROTON[S]','PROTON[I]'},[1 1 -1 1 -1],0,0,1000);
%Phosphate Proton Symport
SymbioModel = addReaction(SymbioModel,'TRSy_Pi',{'Root_Pi[C]','Root_PROTON[C]','Pi[S]','PROTON[S]'},[-1 -2 1 2],0,0,1000);


AATransporters = {'TRSy_GLT','TSyR_L-ASPARTATE','TRSy_VAL','TRSy_LEU','TRSy_ILE','TSyR_L-ALPHA-ALANINE'};
%Proton Pump into the Intermembrane space
SymbioModel = addReaction(SymbioModel,'R_TRI_PROTON',{'Root_ATP[C]','Root_WATER[C]','Root_PROTON[C]','Root_ADP[C]','Root_Pi[C]','PROTON[I]'},[-1 -1 -2 1 1 3],0,0,1000);
%Malate and Succinate to the Rhizobium
SymbioModel = addReaction(SymbioModel,'TRSy_MAL',{'Root_MAL[C]','Root_PROTON[C]','MAL[S]','PROTON[S]'},[-1 -2 1 2],0,0,11.1384);
SymbioModel = addReaction(SymbioModel,'TRSy_SUC',{'Root_SUC[C]','Root_PROTON[C]','SUC[S]','PROTON[S]'},[-1 -2 1 2],0,0,9.4248);

%Hydrogen export and Nitrogen Uptake
SymbioModel = addReaction(SymbioModel,'TSyE_HYDROGEN-MOLECULE',{'HYDROGEN-MOLECULE[S]'},[-1],0,0,1000);
SymbioModel = addReaction(SymbioModel,'TESy_NITROGEN-MOLECULE',{'NITROGEN-MOLECULE[S]'},[1],0,0,1000);
%PHB Storage
SymbioModel = addReaction(SymbioModel,'TSyE_Poly-Hydroxybutyrate',{'Poly-Hydroxybutyrate[S]'},[-1],0,0,1000);
SymbioModel.ub(OxygenPos) = 16.8;
%Turn off the Nitrate and ammonia uptakes
NitrogenUptake = find(ismember(SymbioModel.rxns,{'Root_TEC_NITRATE','Root_TEC_AMMONIUM'}));
SymbioModel.lb(NitrogenUptake) = 0;
SymbioModel.ub(NitrogenUptake) = 0;
%% Do some tests
%Determine the AMMONIA and NITRATE uptake:
%clear all
%load SymbiontConnected
changeCobraSolver('ibm_cplex')
Ammoniumuptake = find(ismember(SymbioModel.rxns,'Root_TEC_AMMONIUM'));
NitrateUptake = find(ismember(SymbioModel.rxns,'Root_TEC_NITRATE'));
%Turn them off and see, whether the model can still grow
SymbioModel.ub([Ammoniumuptake NitrateUptake]) = 0;
SymbioModel.lb([Ammoniumuptake NitrateUptake]) = 0;

ParaModel = addSumFluxConstraint(SymbioModel,{'TRSy_MAL','TRSy_SUC'},[1 1],0,11.1384,'DicarboxExchange');
ParaModel = addSumFluxConstraint(ParaModel,AATransporters,ones(numel(AATransporters),1),0,38.20,'AminoAcidExchange');

solWOMaintenance = parsimOpt(ParaModel);
fprintf('The maximal growth rate without Maintenance is %f\n', solWOMaintenance.objval * 24);
ATPMaintenance = find(ismember(ParaModel.rxns,'RXN-11135_S'));
ParaModel.lb(ATPMaintenance) = 63.47;
sol = parsimOpt(ParaModel);
fprintf('The maximal growth rate with Maintenance is %f\n', sol.objval * 24);
%% Adjust the amount of Asparagine used
% The growth medium contained excessive amounts of nitrogen, thus the
% asparagine fraction (as nitrogen storage) was extraordinarily high, which
% is not to be expected in a symbiotic system.
disp('Adjusting Alanine Amount in the Symbiotic Model');
%clear all
%load ParametrisedModel;
ModelAdjustedForAsn = ParaModel;
RootBiomass = find(ismember(ModelAdjustedForAsn.rxns,'Root_BiomassRoot'));
LeaveBiomass = find(ismember(ModelAdjustedForAsn.rxns,'Leave_BiomassShoot'));
LeaveBiomassWithoutStarch = find(ismember(ModelAdjustedForAsn.rxns,'Leave_BiomassShootWithOutStarch'));
AsparagineLeave = find(ismember(ModelAdjustedForAsn.mets,{'Leave_ASN[C]'}));
AsparagineRoot = find(ismember(ModelAdjustedForAsn.mets,{'Root_ASN[C]'}));
ASNMolWeight = 132.12;

%First do it for the Asparagine without starch
ASNLeaveWOStarchAmounts = ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
AsnLeaveWOStarchWeights = abs(ASNMolWeight/1e6 * ASNLeaveWOStarchAmounts);
AsnLeaveWOStarchWeightChange = AsnLeaveWOStarchWeights - AsnLeaveWOStarchWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch) = 0.05 * ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
ModelAdjustedForAsn.S(:,LeaveBiomassWithoutStarch) = ModelAdjustedForAsn.S(:,LeaveBiomassWithoutStarch) / (1-AsnLeaveWOStarchWeightChange);

% And repeat it for that with starch
ASNLeaveAmounts = ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
AsnLeaveWeights = abs(ASNMolWeight/1e6 * ASNLeaveAmounts);
AsnLeaveWeightChange = AsnLeaveWeights - AsnLeaveWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomass) = 0.05 * ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomass);
ModelAdjustedForAsn.S(:,LeaveBiomass) = ModelAdjustedForAsn.S(:,LeaveBiomass) / (1-AsnLeaveWeightChange);

%And repeat the process for the root
ASNRootAmounts = ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass);
AsnRootWeights = abs(ASNMolWeight/1e6 * ASNRootAmounts);
AsnRootWeightChange = AsnRootWeights - AsnRootWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass) = 0.05 * ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass);
ModelAdjustedForAsn.S(:,RootBiomass) = ModelAdjustedForAsn.S(:,RootBiomass) / (1-AsnRootWeightChange);

SymbioticModel_Adjusted_For_Asn = ModelAdjustedForAsn;
%save('ModelWithAdjustedAsn','CombinedModel','ModelAdjustedForAsn','RootBiomass','ATPMaintenance')

%% And Perform it also for the Combined Model without Symbiont
%clear all
%load ParametrisedModel;
disp('Adjusting Alanine Amount in the Tissue Model');
ModelAdjustedForAsn = Orig_Combined_Model;
RootBiomass = find(ismember(ModelAdjustedForAsn.rxns,'Root_BiomassRoot'));
LeaveBiomass = find(ismember(ModelAdjustedForAsn.rxns,'Leave_BiomassShoot'));
LeaveBiomassWithoutStarch = find(ismember(ModelAdjustedForAsn.rxns,'Leave_BiomassShootWithOutStarch'));
AsparagineLeave = find(ismember(ModelAdjustedForAsn.mets,{'Leave_ASN[C]'}));
AsparagineRoot = find(ismember(ModelAdjustedForAsn.mets,{'Root_ASN[C]'}));
ASNMolWeight = 132.12;

%First do it for the Asparagine without starch
ASNLeaveWOStarchAmounts = ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
AsnLeaveWOStarchWeights = abs(ASNMolWeight/1e6 * ASNLeaveWOStarchAmounts);
AsnLeaveWOStarchWeightChange = AsnLeaveWOStarchWeights - AsnLeaveWOStarchWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch) = 0.05 * ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
ModelAdjustedForAsn.S(:,LeaveBiomassWithoutStarch) = ModelAdjustedForAsn.S(:,LeaveBiomassWithoutStarch) / (1-AsnLeaveWOStarchWeightChange);

% And repeat it for that with starch
ASNLeaveAmounts = ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
AsnLeaveWeights = abs(ASNMolWeight/1e6 * ASNLeaveAmounts);
AsnLeaveWeightChange = AsnLeaveWeights - AsnLeaveWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomass) = 0.05 * ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomass);
ModelAdjustedForAsn.S(:,LeaveBiomass) = ModelAdjustedForAsn.S(:,LeaveBiomass) / (1-AsnLeaveWeightChange);

%And repeat the process for the root
ASNRootAmounts = ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass);
AsnRootWeights = abs(ASNMolWeight/1e6 * ASNRootAmounts);
AsnRootWeightChange = AsnRootWeights - AsnRootWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass) = 0.05 * ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass);
ModelAdjustedForAsn.S(:,RootBiomass) = ModelAdjustedForAsn.S(:,RootBiomass) / (1-AsnRootWeightChange);

AlaAdjustedTissueModel = ModelAdjustedForAsn;

%save('NonSymbioticModelWithAdjustedAsn','CombinedModel')
%% Do calculations delAdjustesd
disp('Setting Final Parameters')
ATPMaintenanceCost = 63.47;
Proton_i = find(ismember(SymbioticModel_Adjusted_For_Asn.mets,'PROTON[I]'));
Proton_s = find(ismember(SymbioticModel_Adjusted_For_Asn.mets,'PROTON[S]'));
ATPSYN = find(ismember(SymbioticModel_Adjusted_For_Asn.rxns,'ATPSYN-RXN_S'));

AsnSol = parsimOpt(SymbioticModel_Adjusted_For_Asn);
fprintf('The maximal growth rate with Adjusted Asparagine is %f\n', AsnSol.objval * 24);
SymbioticModel_Adjusted_For_Asn.lb(ATPMaintenance) = SymbioticModel_Adjusted_For_Asn.lb(ATPMaintenance) *0.5;
AsnSol = parsimOpt(SymbioticModel_Adjusted_For_Asn);
fprintf('The maximal growth rate with lower Maintenance is %f\n', AsnSol.objval * 24);

%Try to alter the ATPSYN reaction
SymbioticModel_Adjusted_For_Asn.lb(ATPMaintenance) = 0;
AsnSol = parsimOpt(SymbioticModel_Adjusted_For_Asn);
fprintf('The maximal growth rate without Maintenance is %f\n', AsnSol.objval * 24);
%restore the maintenance and adjust the Citrate lyase
SymbioticModel_Adjusted_For_Asn.lb(ATPMaintenance) = ATPMaintenanceCost;
CITLyase = find(ismember(SymbioticModel_Adjusted_For_Asn.rxns,'CITLY-RXN_S'));
SymbioticModel_Adjusted_For_Asn.lb(CITLyase) = 0;
AsnSol = parsimOpt(SymbioticModel_Adjusted_For_Asn);
fprintf('The maximal growth rate with Maintenance and irrev Citrate Lyase is %f\n', AsnSol.objval * 24);
SymbioticModel = SymbioticModel_Adjusted_For_Asn;

cd(origDir);
rmpath([scriptPath filesep 'Utilities']);