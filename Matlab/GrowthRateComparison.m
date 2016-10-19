function GrowthRateComparison(SymbioticModel,CombinedModel)
% Make a comparison Scan for different Availabilities of Nitrate in the
%symbiotic and non Symbiotic models, with a maximal carbon amount
% either use the outputs from BuildSymbioticModel as inputs or run the
% function wihout arguments.
%Using Tisuemodel gives a skewed result, as the biomasses in both models
%are different due to the adjustment in asparagine in the Symbiotic
%condition.

%restricted by the Overnight starch amount, which is 0.0125*0.66g/g  or 
% 46.32 umol/g/h
% 5.787037037
%and its night 8/24 hours per day
umol_Starch_per_h = 0.0125 * 2/3 / 180 *1e6 * 8/24;
%Average daily Light uptake is assumed to be 1000 Photons/gDW/h
Lightuptake = 1000 * 16/24; 
%So this is our average growth.

%Create the models from Scratch

[SymbioticModel,CombinedModel] = BuildSymbioticModel();


%First, get the Light uptake and close it
%load NonSymbioticModelWithAdjustedAsn
%load ModelWithAsnAndDeactCitLy
%Set up the non symbiotic model
%Shut off Light
LightImport = find(ismember(CombinedModel.rxns,'Leave_TEH_Light'));
CombinedModel.ub(LightImport) = Lightuptake;
%Set starch to the amount available over night
StarchImport = find(ismember(CombinedModel.rxns,'Leave_TEH_Starch'));
CombinedModel.ub(StarchImport) = umol_Starch_per_h;
CombinedModel.c(:) = 0;
BiomassNonSym = find(ismember(CombinedModel.rxns,'Biomass'));
LeaveBiomass= find(ismember(CombinedModel.rxns,'Leave_BiomassShoot'));
CombinedModel.ub(LeaveBiomass) = 0;
CombinedModel.c(BiomassNonSym) = 1;
NitrateUptakeNonSym = find(ismember(CombinedModel.rxns,'Root_TEC_NITRATE'));
NitrateUptakeSym = find(ismember(SymbioticModel.rxns,'Root_TEC_NITRATE'));
CombinedModel.ub(NitrateUptakeNonSym) = 0;
%Also impose the growth Maintenance.
%This derives from the numbers for Trifolum
%There: 15 mg Glucose/gDW/day are consumed.
%1 mol of Starch allows the regeneration of 32 mol of ATP
%15mg/24h = 0.6250  mg/gDW h
%0.6250 mg/ gDW h / 180 g/mol  * 32 mol ATP/mol=  0.1111 mmol ATP/ gDW h
% or in umol = 111.1 umol/gDW/h
%We will split this to Root (1/3) and shoot  (2/3).
MaintenanceRequirement = 0.6250 / 180 * 32 * 1000;

ATPMaintenanceShoot = find(ismember(CombinedModel.rxns,'Leave_ATPSYN-RXN_C'));
ATPMaintenanceRoot = find(ismember(CombinedModel.rxns,'Root_ATPSYN-RXN_C'));
CombinedModel.lb(ATPMaintenanceShoot) = MaintenanceRequirement* 2/3;
CombinedModel.lb(ATPMaintenanceRoot) = MaintenanceRequirement* 1/3;

%Set up the symbiotic model
LightImport = find(ismember(SymbioticModel.rxns,'Leave_TEH_Light'));
SymbioticModel.ub(LightImport) = Lightuptake;
StarchImport = find(ismember(SymbioticModel.rxns,'Leave_TEH_Starch'));
SymbioticModel.ub(StarchImport) = umol_Starch_per_h;
SymbioticModel.c(:) = 0;
BiomassSym = find(ismember(SymbioticModel.rxns,'Biomass'));
LeaveBiomass= find(ismember(SymbioticModel.rxns,'Leave_BiomassShoot'));
SymbioticModel.ub(LeaveBiomass) = 0;
SymbioticModel.c(BiomassSym) = 1;
ATPMaintenanceShoot = find(ismember(SymbioticModel.rxns,'Leave_ATPSYN-RXN_C'));
ATPMaintenanceRoot = find(ismember(SymbioticModel.rxns,'Root_ATPSYN-RXN_C'));
SymbioticModel.lb(ATPMaintenanceShoot) = MaintenanceRequirement* 2/3;
SymbioticModel.lb(ATPMaintenanceRoot) = MaintenanceRequirement* 1/3;



solSym = parsimOpt(SymbioticModel);
solNonSym = parsimOpt(CombinedModel);

disp(['The Maximal growth rate for non symbiotic conditions is ' num2str(24* solNonSym.objval)]);
disp(['The Maximal growth rate for symbiotic conditions is ' num2str(24* solSym.objval)]);

Ammoniumuptake = find(ismember(CombinedModel.rxns,'Root_TEC_AMMONIUM'));
AmmoniumuptakeSym = find(ismember(SymbioticModel.rxns,'Root_TEC_AMMONIUM'));
CombinedModel.ub(NitrateUptakeNonSym) = 0;
CombinedModel.ub(Ammoniumuptake) = 0;
SymbioticModel.ub(AmmoniumuptakeSym) = 0;
SymbioticModel.ub(NitrateUptakeNonSym) = 0;
nonsymscans = ScanUpperBound(CombinedModel,'Root_TEC_AMMONIUM',0,25,30);
symscans = ScanUpperBound(SymbioticModel,'Root_TEC_AMMONIUM',0,25,30);
Vals = linspace(0,25,30);

%% And plot the Figure
f = figure
p = plot(Vals,[nonsymscans(BiomassNonSym,:)*24;symscans(BiomassSym,:)*24],'LineWidth',3); % ,ATPSynthase_Root,ATPSynthase_Shoot,
LineStyles = {'-','-.'};
for i=1:min(numel(LineStyles))
p(i).LineStyle = LineStyles{i};
end
caxes = gca;
xlabel('Available ammonium ($\frac{\mu mol}{g \cdot h}$)','Interpreter','latex','FontSize',20);
ylabel('Maximal biomass yield ($\frac{mg}{g \cdot h}$)','Interpreter','latex','FontSize',20);
leg = legend({'Non symbiotic system','Symbiotic system'});
leg.Interpreter = 'latex';
leg.Position = [0.5 0.45 0.3 0.2];
set(caxes,'FontSize',18);
set(caxes,'Position',[0.2 0.19 0.77 0.74],'LineWidth',2,'TickLabelInterpreter','latex','XLim',[0 25]);

