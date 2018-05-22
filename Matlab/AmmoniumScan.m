function TissueModel = AmmoniumScan(CombinedModel)
%% Scan Ammonia
if nargin == 0
    CombinedModel = BuildTissueModel();
end

%adjust pathes
scriptPath = fileparts(which(mfilename));
origDir = cd(scriptPath);
addpath([scriptPath filesep 'Utilities']);

TissueModel = CombinedModel;

%% Get the positions of the various Reactions and exchangers
WaterExchange = find(ismember(CombinedModel.rxns,{'Root_TEC_WATER','Leave_TCE_WATER','TSR_WATER','TRS_WATER'}));
CombinedModel.lb(WaterExchange(1)) = min(CombinedModel.lb);
CO2Exchange = find(ismember(CombinedModel.rxns,{'Root_TCE_CARBON-DIOXIDE','Leave_TCE_CARBON-DIOXIDE'}));
HCO3Conversion = find(ismember(CombinedModel.rxns,{'Root_RXN0-5224_C','Leave_RXN0-5224_C'}));
LightImport = find(ismember(CombinedModel.rxns,'Leave_TEH_Light'));
StarchImport = find(ismember(CombinedModel.rxns,'Leave_TEH_Starch'));
RootBiomass = find(ismember(CombinedModel.rxns,'Root_BiomassRoot'));
Biomass = find(ismember(CombinedModel.rxns,'Biomass'));
LeaveBiomassWithOutStarch = find(ismember(CombinedModel.rxns,'Leave_BiomassShootWithOutStarch'));
LeaveBiomass = find(ismember(CombinedModel.rxns,'Leave_BiomassShoot'));
NitrateUptake = find(ismember(CombinedModel.rxns,'Root_TEC_NITRATE'));
Ammoniumuptake = find(ismember(CombinedModel.rxns,'Root_TEC_AMMONIUM'));
ProtonExchange = find(ismember(CombinedModel.rxns,'Root_TCE_PROTON'));
GlycolysisStartRoot = find(ismember(CombinedModel.rxns,'Root_F16ALDOLASE-RXN_C'));
GlycolysisStartShoot = find(ismember(CombinedModel.rxns,'Leave_F16ALDOLASE-RXN_C'));
GlycolysisEndRoot = find(ismember(CombinedModel.rxns,'Root_GAPOXNPHOSPHN-RXN_C'));
GlycolysisEndShoot = find(ismember(CombinedModel.rxns,'Leave_GAPOXNPHOSPHN-RXN_C'));
TCA_root = find(ismember(CombinedModel.rxns,'Root_SUCCCOASYN-RXN_M'));
TCA_shoot = find(ismember(CombinedModel.rxns,'Leave_SUCCCOASYN-RXN_M'));
PPP_root = find(ismember(CombinedModel.rxns,'Root_GLU6PDEHYDROG-RXN_H'));
PPP_shoot = find(ismember(CombinedModel.rxns,'Leave_GLU6PDEHYDROG-RXN_H'));
HCO3Exchange = find(ismember(CombinedModel.rxns,'Root_TCE_HCO3'));
LeaveHCO3Exchange = find(ismember(CombinedModel.rxns,'Leave_TCE_HCO3'));
LeaveGlucoseStorage = find(ismember(CombinedModel.rxns,'Leave_TCE_GLC'));
RootValineStorage = find(ismember(CombinedModel.rxns,'Root_TCE_VAL'));
LeaveSuccinateStorage = find(ismember(CombinedModel.rxns,'Leave_TCE_SUC'));
ATPSynthase_Shoot = find(ismember(CombinedModel.rxns,'Leave_ATPSYN-RXN_M'));
ATPSynthase_Root= find(ismember(CombinedModel.rxns,'Root_ATPSYN-RXN_M'));

%Impose Night Coniditions and don't allow Glucose storage:
CombinedModel.lb([LeaveHCO3Exchange,LightImport,LeaveBiomass]) = 0;
CombinedModel.ub([LeaveHCO3Exchange,LightImport,LeaveBiomass]) = 0;

%Also, don'T allow internal storage
CombinedModel.lb([LeaveSuccinateStorage]) = 0;
CombinedModel.ub([LeaveSuccinateStorage]) = 0;


Glycolysis = {'PGLUCISOM-RXN',{'6PFRUCTPHOS-RXN','2.7.1.90-RXN'},'F16ALDOLASE-RXN','GAPOXNPHOSPHN-RXN','PHOSGLYPHOS-RXN','3PGAREARR-RXN','2PGADEHYDRAT-RXN','PEPDEPHOS-RXN'};
GlycRootPositions = {};
GlycShootPositions = {};
for i = 1:numel(Glycolysis)
    rootpos = [];
    shootpos = [];
    if iscell(Glycolysis{i})
        current = Glycolysis{i};
        for j = 1:numel(current)
            rootpos = [rootpos ; find(~cellfun(@isempty , strfind(CombinedModel.rxns,current{j})) & ~cellfun(@isempty , strfind(CombinedModel.rxns,'Root')))];
            shootpos = [shootpos ; find(~cellfun(@isempty , strfind(CombinedModel.rxns,current{j})) & ~cellfun(@isempty , strfind(CombinedModel.rxns,'Leave')))];
        end
    else
        rootpos = find(~cellfun(@isempty , strfind(CombinedModel.rxns,Glycolysis{i}))& ~cellfun(@isempty , strfind(CombinedModel.rxns,'Root')));
        shootpos = find(~cellfun(@isempty , strfind(CombinedModel.rxns,Glycolysis{i}))& ~cellfun(@isempty , strfind(CombinedModel.rxns,'Leave')));
    end
    GlycRootPositions{i} = rootpos;
    GlycShootPositions{i} = shootpos;
end


%Fix Biomass production to some "reasonable" amounts (i.e. the normal growth rate of 0.1 g/gDW)
%since our Biomass is adjusted to be 1g/Flux (given the flux is in umol/gDW h)
%we can simply set it to 0.1 / 24 and assume a constant growth during day
%and night
GrowthPerHour = 0.1/24;
CombinedModel.lb(Biomass) = GrowthPerHour;
CombinedModel.ub(Biomass) = GrowthPerHour;

%Also impose the growth Maintenance.
%This derives from the numbers for Trifolum
%There: 15 mg Glucose/gDW/day are consumed.
%1 mol of Starch allows the regeneration of 32 mol of ATP
%15mg/24h = 0.6250  mg/gDW h
%0.6250 mg/ gDW h / 180 g/mol  * 32 mol ATP/mol=  0.1111 mmol ATP/ gDW h
% or in umol = 111.1 umol/gDW/h
%We will split this to Root (1/3) and shoot  (2/3).
ATPMaintenanceShoot = find(ismember(CombinedModel.rxns,'Leave_ATPSYN-RXN_C'));
ATPMaintenanceRoot = find(ismember(CombinedModel.rxns,'Root_ATPSYN-RXN_C'));
CombinedModel.lb(ATPMaintenanceShoot) = 111.1* 2/3;
CombinedModel.lb(ATPMaintenanceRoot) = 111.1* 1/3;

%Obtain the minimal uptake rate of Ammonium, when no Nitrate is available
%(i.e. the maximal Ammonium use)
CombinedModel.c(:) = 0;
CombinedModel.c(Ammoniumuptake) = -1;
CombinedModel.ub(NitrateUptake) = 0;
symsol = parsimOpt(CombinedModel);
AmmoniumMax = symsol.x(Ammoniumuptake);


minweight = 0;
maxweight = AmmoniumMax;
CombinedModel.c(:) = 1;
%Allow free exchange of Water, export and uptake of CO2 and Conversion of
%H+ + HCO3 to CO2 + H2O
CombinedModel.c(WaterExchange) = 0;
CombinedModel.c(CO2Exchange) = 0;
CombinedModel.c(HCO3Conversion) = 0;
Steps = 20;
%%
CombinedModel.ub(NitrateUptake) = max(CombinedModel.ub);

%Perform the Scans
scans = ScanDualReactionQP(CombinedModel,CombinedModel.rxns(Ammoniumuptake),CombinedModel.rxns(NitrateUptake),minweight,maxweight,Steps);

%Our TCA reaction is a reaction that normally runs in reverse, so just flip
%the sign in the results.
scans([TCA_shoot,TCA_root],:) = -scans([TCA_shoot,TCA_root],:);
scans(end+1,:) = scans(ProtonExchange,:) - scans(HCO3Exchange,:);
ProtonExchangePos = size(scans,1);

%And combine the root and shoot TCA.
scans(end+1,:) = scans(TCA_shoot,:) + scans(TCA_root,:);
TCAPos = size(scans,1);
for i = 1:numel(Glycolysis)
    scans(end+1,:) = sum(scans(GlycRootPositions{i},:));
end
for i = 1:numel(Glycolysis)
    scans(end+1,:) = sum(scans(GlycShootPositions{i},:));
end
Weights = linspace(minweight,maxweight,Steps);
%%
close all
%Get some positions, that we might want to plot.
FluxPos = [StarchImport, NitrateUptake,PPP_root,PPP_shoot,HCO3Exchange,ProtonExchange,GlycolysisEndRoot,GlycolysisEndShoot,GlycolysisStartRoot,GlycolysisStartShoot,ATPSynthase_Root,ATPSynthase_Shoot,TCA_shoot,TCA_root,ProtonExchangePos, TCAPos];
FluxName = {'Starch consumption','Nitrate uptake','PPP in root','PPP in shoot','HCO$_3^-$export','Proton export','Glycolysis End Root','Glycolysis End Shoot','Glycolysis Start Root','Glycolysis Start Shoot','ATPSynthase in root', 'ATPSynthase in shoot','TCA in shoot','TCA in root','Charge Balance','TCA'};
Starch = 1;
Nitrate = 2;
Root_PPP = 3;
Shoot_PPP = 4;
HCO3 = 5;
Protons = 6;
GlycolysisPayRoot = 7;
GlycolysisPayShoot = 8;
GlycolysisPrepRoot = 9;
GlycolysisPrepShoot = 10;
ATPSynRoot = 11;
ATPSynShoot = 12;
TCAShoot = 13;
TCARoot = 14;
ChargeBalance = 15;
TCA = 16;
plotteditems = [Starch,Nitrate,Protons,HCO3,TCA, ATPSynShoot];

%Set up the figure and plot the scans.
f = figure;
f.Position = [100 100 1200 600];
p = plot(Weights,scans(FluxPos(plotteditems),:),'LineWidth',3); 
caxes = gca;
paper = gcf;
paper.PaperSize = [5.5 2.5];
paper.PaperPosition = [0 0 5.5 2.5];
set(caxes,'Position',[0.1 0.13 0.57 0.82],'LineWidth',2,'TickLabelInterpreter','latex','XLim',[0,maxweight])
LineStyles = {'-.','-.','--','--',':',':','-','-','-.','-.','--','--',':',':','-','-'};
for i=1:min(numel(LineStyles),numel(plotteditems))
    p(i).LineStyle = LineStyles{i};
end
%Create Legends and axis labels
leg = legend(FluxName(plotteditems));
leg.FontSize = 20;
leg.Position = [0.70 0.33 0.25  0.4];
leg.Interpreter = 'latex';
yl = ylabel('Fluxes in ($\mu mol~g^{-1} h^{-1}$)','Interpreter','latex','FontSize',20);
xl = xlabel('Ammonium uptake flux in ($\mu mol~g^{-1} h^{-1}$)','Interpreter','latex','FontSize',20);

cd(origDir);
rmpath([scriptPath filesep 'Utilities']);