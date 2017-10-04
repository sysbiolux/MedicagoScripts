function CombinedModel = BuildTissueModel(Medicago)
%Build the Tissue Model based on the XML File of the Medicago Model.
if nargin == 1    
    OrigMedicago = Medicago;
else
    OrigMedicago = importMedicago();
end
    
OrigMedicago.ub(ismember(OrigMedicago.rxns,'TCE_GLC')) = 0;
OrigMedicago.lb(ismember(OrigMedicago.rxns,'TCE_GLC')) = 0;
RootPresence = importdata(['Data' filesep 'RootPresenceCalls.csv']);
RootPresence = RootPresence.data;
StemPresence = importdata(['Data' filesep 'StemPresenceCalls.csv']);
StemPresence = StemPresence.data;
LeavePresence = importdata(['Data' filesep 'LeavePresenceCalls.csv']);
LeavePresence = LeavePresence.data;
%This import is matched to the respective rows in the Presence Calls.
load(['Data' filesep 'GeneIDForPresenceCall.mat'])

%Benedito2009:
%2+ Present Calls  -> Pesent -> is 1 in the data
%1 Present Call  -> Unknown -> is 0 in the used data
%0 Present Calls -> Absent -> is -1 in the used data


%Steps for simulation:
%Determine a Tissue specific network using fastcore and all present reactions as core.
%In addition, the following reactions are set as part of the core:
%1. Nitrate uptake
%2. Ammonium uptake
%3. Biomass (with Starch for leaves, without starch for root).
%Create a splitted LP. 
%Modify the weight on the ammonium uptake


%% Create a simultaneous model
% Create Root and Shoot
%clear all
%load BaseData2;
Medicago = OrigMedicago

%Remove ATPase and DEHOG
ATPase = find(ismember(Medicago.rxns,'ATPase'));
DEHOG = find(ismember(Medicago.rxns,'DEHOG'));


Medicago = removeRxns(Medicago,Medicago.rxns([ATPase,DEHOG]));
Root = Medicago;
Leave = Medicago;



%Remove Unavailable uptakes
%These are light and all Carbon sources for the Root.
%Also Turn off the Glucose Export
disp('Removing Root Reactions')
RootRemovedReacs = find(ismember(Root.rxns,{'TEC_GLC','TEC_SUCROSE','TEH_Light','TEH_Starch'}));
Root = removeRxns(Root,Root.rxns(RootRemovedReacs));
Root = removeRxns(Root,'TCE_GLC');
%And any other uptake for the shoot
% Only retain CO2 , O2, Light (Day) and Starch (Night)
%We retain the remaining exporters as those compounds could be stored.
%The only exceptions are: PROTONS and Glucose
disp('Removing Importers from Shoot')
ImportersAndExportersRemovedInLeave = find(ismember(Leave.rxns,{'TEC_AMMONIUM',    'TEC_CO+2',    'TEC_CPD-3',    'TEC_FE+2',    'TEC_GLC',    'TEC_MG+2',    'TEC_NITRATE',    'TEC_Pi',    'TEC_SUCROSE',    'TEC_SULFATE',    'TEC_WATER', 'TCE_VAL'}));
Leave = removeRxns(Leave,Leave.rxns(ImportersAndExportersRemovedInLeave)) ;
Leave = removeRxns(Leave,'TCE_PROTON');

Root.rxns = strcat('Root_', Root.rxns);
Root.mets = strcat('Root_', Root.mets);

%Root:
Leave.rxns = strcat('Leave_', Leave.rxns);
Leave.mets = strcat('Leave_', Leave.mets);

disp('Combining Models');
CombinedModel = combineModels(Leave,Root);


ShootUptake = {'GLN','L-ASPARTATE','GLT','NITRATE','L-ALPHA-ALANINE','PROTON','WATER','SULFATE',...
              'MG+2','FE+2','Pi','CO+2','CPD-3'};
%Gln, Asp, Glt, Ala, Nitrate, Water, Protons, SO4, Mg, Fe,
%Ammonia,Molybedenium, Cobolt, 
%Phophate
RootDelivery = {'ASN','L-ASPARTATE','L-ALPHA-ALANINE','GLN','GLT','SER',...
                'PRO','GLY','THR','VAL','ILE','LEU','LYS','ARG','HIS','WATER',...
                'PROTON','SUCROSE'};
%Towards Root:
%Asn, Asp, Ala, Gln, Glt, Ser, Pro, Gly, Thr, Val, Ile, Leu, Lys, Arg, His,
%Water, Protons
disp('Adding Transporters')
%Ammonium as cotransport!
NoTransport = CombinedModel.rxns;
CombinedModel = addReaction(CombinedModel,'TRS_AMMONIUM',{'Root_AMMONIA[C]','Root_PROTON[C]','Leave_PROTON[C]','Leave_AMMONIA[C]'},[-1 -1 1 1],0,0,1000);
for i=1:numel(ShootUptake)
    CombinedModel = addReaction(CombinedModel,['TRS_' ShootUptake{i}],{['Root_' ShootUptake{i} '[C]'],['Leave_' ShootUptake{i} '[C]'] }, [-1 1], 0,0,1000);
end
for i=1:numel(RootDelivery)
    CombinedModel = addReaction(CombinedModel,['TSR_' RootDelivery{i}],{['Leave_' RootDelivery{i} '[C]'],['Root_' RootDelivery{i} '[C]'] }, [-1 1], 0,0,1000);
end

Transporters = setdiff(CombinedModel.rxns,NoTransport);
% Determine 4 Models (One Day Model, One Night model)
% Day Model: Allow Light, No Starch, BiomassWithStarch for Leave, RootBiomass for Root
% Night Model: No Light, Allow Starch, Biomass for Leave, RootBiomass for Root
% For Day and night use Ammonia and Nitrate for growth.
LightImport = find(ismember(CombinedModel.rxns,'Leave_TEH_Light'));
StarchImport = find(ismember(CombinedModel.rxns,'Leave_TEH_Starch'));
RootBiomass = find(ismember(CombinedModel.rxns,'Root_BiomassRoot'));
LeaveBiomassWithStarch = find(ismember(CombinedModel.rxns,'Leave_BiomassShoot'));
LeaveBiomass = find(ismember(CombinedModel.rxns,'Leave_BiomassShootWithOutStarch'));
NitrateUptake = find(ismember(CombinedModel.rxns,'Root_TEC_NITRATE'));
Ammoniumuptake = find(ismember(CombinedModel.rxns,'Root_TEC_AMMONIUM'));

%Determine Gene Positions
[A,B] = ismember(GeneID,Medicago.genes);

%Determine the Leave Present Reactions
LeaveActivities = zeros(1,numel(Medicago.genes));
LeaveActivities(B(B~=0)) = LeavePresence(A);
LeaveOn = Medicago.genes(LeaveActivities==1);
StemActivities = zeros(1,numel(Medicago.genes));
StemActivities(B(B~=0)) = StemPresence(A);
StemOn = Medicago.genes(StemActivities==1);


LeaveReacs = findRxnsActiveWithGenes(Medicago,LeaveOn);
StemReacs = findRxnsActiveWithGenes(Medicago,StemOn);
ShootReacs = union(LeaveReacs,StemReacs);

LeaveReacs = strcat('Leave_',ShootReacs);
%Determine The Root Present Reactions
RootActivities = zeros(1,numel(Medicago.genes));
RootActivities(B(B~=0)) = RootPresence(A);
RootOn = Medicago.genes(RootActivities==1);
RootReacs = findRxnsActiveWithGenes(Medicago,RootOn);
RootReacs = strcat('Root_',RootReacs);
CoreReacs = union(RootReacs,LeaveReacs);
DayAmmoniaReacs = union(CoreReacs,CombinedModel.rxns([Ammoniumuptake,LightImport,LeaveBiomassWithStarch,RootBiomass]));
DayNitrateReacs = union(CoreReacs,CombinedModel.rxns([NitrateUptake,LightImport,LeaveBiomassWithStarch,RootBiomass]));
NightAmmoniaReacs = union(CoreReacs,CombinedModel.rxns([Ammoniumuptake,StarchImport,LeaveBiomass,RootBiomass]));
NightNitrateReacs = union(CoreReacs,CombinedModel.rxns([NitrateUptake,StarchImport,LeaveBiomass,RootBiomass]));

disp('Making Model Consistent')
ConsistentCombined = fastcc(CombinedModel,1e-4);

ConsistentModel = removeRxns(CombinedModel,setdiff(CombinedModel.rxns,CombinedModel.rxns(ConsistentCombined)));

%save('ConsistentFromCombined_From_XML','CombinedModel','ConsistentModel','DayAmmoniaReacs','DayNitrateReacs','NightAmmoniaReacs','NightNitrateReacs','Transporters');

%%
%clear all
%load ConsistentFromCombined
%load BaseData2

disp('Generating Models for different Conditions')
LightImport = find(ismember(ConsistentModel.rxns,'Leave_TEH_Light'));
StarchImport = find(ismember(ConsistentModel.rxns,'Leave_TEH_Starch'));
RootBiomass = find(ismember(ConsistentModel.rxns,'Root_BiomassRoot'));
LeaveBiomassWithStarch = find(ismember(ConsistentModel.rxns,'Leave_BiomassShoot'));
LeaveBiomass = find(ismember(ConsistentModel.rxns,'Leave_BiomassShootWithOutStarch'));
NitrateUptake = find(ismember(ConsistentModel.rxns,'Root_TEC_NITRATE'));
Ammoniumuptake = find(ismember(ConsistentModel.rxns,'Root_TEC_AMMONIUM'));

DayAmmoniumModel = ConsistentModel;
DayAmmoniumModel.ub([StarchImport,LeaveBiomass,NitrateUptake]) = 0;

DayNitrateModel = ConsistentModel;
DayNitrateModel.ub([StarchImport,LeaveBiomass,Ammoniumuptake]) = 0;

NightAmmoniumModel = ConsistentModel;
NightAmmoniumModel.ub([LightImport,LeaveBiomassWithStarch,NitrateUptake]) = 0;

NightNitrateModel = ConsistentModel;
NightNitrateModel.ub([LightImport,LeaveBiomassWithStarch,Ammoniumuptake]) = 0;


%% Make all those models consistent
disp('Making Day Ammonium model Consistent')
A = fastcc(DayAmmoniumModel,1e-4);
DayAmmoniumModel = removeRxns(DayAmmoniumModel,setdiff(DayAmmoniumModel.rxns,DayAmmoniumModel.rxns(A)));

disp('Making Day Nitrate model Consistent')
A = fastcc(DayNitrateModel,1e-4);
DayNitrateModel = removeRxns(DayNitrateModel,setdiff(DayNitrateModel.rxns,DayNitrateModel.rxns(A)));

disp('Making Night Ammonium model Consistent')
A = fastcc(NightNitrateModel,1e-4);
NightNitrateModel = removeRxns(NightNitrateModel,setdiff(NightNitrateModel.rxns,NightNitrateModel.rxns(A)));

disp('Making Night nitrate model Consistent')
A = fastcc(NightAmmoniumModel,1e-4);
NightAmmoniumModel = removeRxns(NightAmmoniumModel,setdiff(NightAmmoniumModel.rxns,NightAmmoniumModel.rxns(A)));

%% And Launch Fastcore for the Models
%
CoreDayAmmonium = find(ismember(DayAmmoniumModel.rxns,DayAmmoniaReacs));
NoPenaltyDayAmmonium = find(ismember(DayAmmoniumModel.rxns,Transporters));
CoreDayNitrate = find(ismember(DayNitrateModel.rxns,DayNitrateReacs));
NoPenaltyDayNitrate = find(ismember(DayNitrateModel.rxns,Transporters));
CoreNightAmmonium = find(ismember(NightAmmoniumModel.rxns,NightAmmoniaReacs));
NoPenaltyNightAmmonium = find(ismember(NightAmmoniumModel.rxns,Transporters));
CoreNightNitrate = find(ismember(NightNitrateModel.rxns,NightNitrateReacs));
NoPenaltyNightNitrate = find(ismember(NightNitrateModel.rxns,Transporters));
%%
disp('Using fastcore to restrict Day Ammonium Model')
DayAmmoniumModel = fastcore(DayAmmoniumModel, CoreDayAmmonium, 1e-4);
%DayAmmoniumModel = removeRxns(DayAmmoniumModel,setdiff(DayAmmoniumModel.rxns,DayAmmoniumModel.rxns(A)));

disp('Using fastcore to restrict Day Nitrate Model')
DayNitrateModel = fastcore(DayNitrateModel,CoreDayNitrate,1e-4);
%DayNitrateModel = removeRxns(DayNitrateModel,setdiff(DayNitrateModel.rxns,DayNitrateModel.rxns(A)));

disp('Using fastcore to restrict Night Nitrate Model')
NightNitrateModel = fastcore(NightNitrateModel, CoreNightNitrate, 1e-4);
%NightNitrateModel = removeRxns(NightNitrateModel,setdiff(NightNitrateModel.rxns,NightNitrateModel.rxns(A)));

disp('Using fastcore to restrict Night Ammonium Model')
NightAmmoniumModel = fastcore(NightAmmoniumModel, CoreNightAmmonium, 1e-4);
%NightAmmoniumModel = removeRxns(NightAmmoniumModel,setdiff(NightAmmoniumModel.rxns,NightAmmoniumModel.rxns(A)));

disp('Adding Combined Biomass Reaction')
%Combine All Models
CombinedModelReacs = union(Transporters,union(union(NightAmmoniumModel.rxns,NightNitrateModel.rxns),union(DayAmmoniumModel.rxns,DayNitrateModel.rxns)));
%Keep the Root Proton Exchanger
CombinedModelReacs{end+1} = 'Root_TCE_PROTON';
CombinedModel = removeRxns(CombinedModel,setdiff(CombinedModel.rxns,CombinedModelReacs));
%Now, add a reaction combining Biomass for root and shoot
CombinedModel = addReaction(CombinedModel,'Biomass',{'BiomassRoot','BiomassShoot'} ,[-1 -1],0,0,1000);
ShootBiomassReactions = find(ismember(CombinedModel.rxns,{'Leave_BiomassShoot','Leave_BiomassShootWithOutStarch'}));
RootBiomassReactions = find(ismember(CombinedModel.rxns,{'Root_BiomassRoot'}));
BiomassShootPos = find(ismember(CombinedModel.mets,'BiomassShoot'));
BiomassRootPos = find(ismember(CombinedModel.mets,'BiomassRoot'));
CombinedModel.S(BiomassShootPos,ShootBiomassReactions) = 1;
CombinedModel.S(BiomassRootPos,RootBiomassReactions) = 1;
CombinedModel.S(BiomassShootPos,end) = -2/3;
CombinedModel.S(BiomassRootPos,end) = -1/3;
CombinedModel.c(:) = 0;
CombinedModel.c(end) = 1;
CombinedModel.ub(:) = CombinedModel.ub(:)*1000;
CombinedModel.lb(:) = CombinedModel.lb(:)*1000;

disp('Testing the model')

changeCobraSolver('ibm_cplex');
optimizeCbModel(CombinedModel)
CombinedModel = rmfield(CombinedModel,'modelVersion')
