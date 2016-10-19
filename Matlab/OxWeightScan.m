function OxWeightScan(SymbioticModel)
%Make a comparison Scan for different Availabilities of Nitrate in the
%symbiotic and non Symbiotic models, with a maximal carbon amount
%restricted by the Overnight starch amount, which is 0.0125*0.66g/g  or 
if nargin == 0
    SymbioticModel = BuildSymbioticModel();
end

maxscanval = 20;
scanpoints = 50;
maxOxUptake = 19.76;
SymbioticBiomass = 0.1*0.8/24;



%The normal growth rate of Medicago is at ~0.1 g/gDay
%Symbiotic growth is commonly somewhere between 60 and 90% of this
%So we set it to 0.1/24*0.6, which is in the lower range of those values to allow more flexibility.

LightImport = find(ismember(SymbioticModel.rxns,'Leave_TEH_Light'));
StarchImport = find(ismember(SymbioticModel.rxns,'Leave_TEH_Starch'));

atpprod = find(ismember(SymbioticModel.rxns,{'PEPCARBOX-RXN_S','PYRUVATE-CARBOXYLASE-RXN_S'}));


SymbioticModel.ub(StarchImport) = 0;
SymbioticModel.c(:) = 0;
BiomassSym = find(ismember(SymbioticModel.rxns,'Biomass'));
SymbioticModel.lb(BiomassSym) = SymbioticBiomass;
LeaveBiomassWithoutStarch= find(ismember(SymbioticModel.rxns,'Leave_BiomassShootWithOutStarch'));
SymbioticModel.ub(LeaveBiomassWithoutStarch) = 0;
OxygenPos = find(ismember(SymbioticModel.rxns,'TRSy_OXYGEN-MOLECULE')); 
ATPMaintenance = find(ismember(SymbioticModel.rxns,'Root_ATPSYN-RXN_C'));
%125 umol is about the ammount derived from 15mg glucose /day (5mg/night)
ATPMaintenanceShoot = find(ismember(SymbioticModel.rxns,'Leave_ATPSYN-RXN_C'));
ATPMaintenanceRoot = find(ismember(SymbioticModel.rxns,'Root_ATPSYN-RXN_C'));
SymbioticModel.lb(ATPMaintenanceShoot) = 111.1* 2/3;
SymbioticModel.lb(ATPMaintenanceRoot) = 111.1* 1/3;
%SymbioticModel.ub(ismember(SymbioticModel.rxns,'TSyE_Poly-Hydroxybutyrate')) = 0;
SymbioticModel.lb(ismember(SymbioticModel.rxns,'TSyR_L-ASPARTATE')) = 0;
SymbioticModel.ub(OxygenPos) = maxOxUptake;
SymbioticModel.c(:) = 0;
SymbioticModel.c(OxygenPos) = -1;
symsol = parsimOpt(SymbioticModel);
SymbioticModel.c(:) = 1;
symsol2 = parsimOptQP(SymbioticModel);
optOx = symsol2.x(OxygenPos);

%% Do Scans

symscans = ScanSingleReactionQP(SymbioticModel,'TRSy_OXYGEN-MOLECULE',symsol.x(OxygenPos),maxscanval,scanpoints);



%%
ReacIDs = {'TSyR_AMMONIA';'TSyR_L-ALPHA-ALANINE';'TRSy_VAL';'TRSy_LEU';'TRSy_ILE';'TRSy_OXYGEN-MOLECULE';'TRSy_GLT';'TSyR_L-ASPARTATE';'TSyR_AMMONIUM';'TRSy_MAL';'TRSy_SUC';'TSyE_Poly-Hydroxybutyrate'};
RxnNames = {'Ammonia export';'Alanine export';'Valine uptake';'Leucine uptake';'Isoleucine uptake';'Oxygen uptake';'Glutamate uptake';'Aspartate exchange';'Ammonium export';'Malate uptake';'Succinate uptake';'PHB production'};

AminoAcidUsage =  [findRxnsFromMets(SymbioticModel,{'GLT[S]','2-KETOGLUTARATE[S]'},'containsAll',1) ; findRxnsFromMets(SymbioticModel,{'L-ALPHA-ALANINE[S]'},'containsAll',1) ;  findRxnsFromMets(SymbioticModel,{'VAL[S]'},'containsAll',1)]; 

CarboxylateUsage = [findRxnsFromMets(SymbioticModel,{'OXALACETIC_ACID[S]'},'containsAll',1)];
%ReacIDs = [ReacIDs ; CarboxylateUsage];

%RxnNames = [RxnNames ; CarboxylateUsage];



exchangers  = ismember(SymbioticModel.rxns, ReacIDs);
%use the stds to see if there is any change in the values)
stds = std(symscans')';
relexch = exchangers & abs(stds) > 1e-3;
relexch(OxygenPos) = 0;
[rpres,rpos] = ismember(SymbioticModel.rxns(relexch),ReacIDs);
relreacs = rpos(rpres);
%% Create Figure
%close all
f = figure;
f.Position = [50 50 1200 600];
Vals = linspace(symsol.x(OxygenPos),maxscanval,scanpoints);
%normscans = normalizeScans(symscanns);
p = plot(Vals,symscans(relexch,:),'LineWidth',3); 
LineStyles = {'-','-','-.','-.','-.','--','--',':',':',':','-.','-.','--','--',':',':','-','-'};
LineStyles = [LineStyles, {'-','-.','-.','-.','--','--',':',':',':','-.','-.','--','--',':',':','-','-'}];
LineStyles = [LineStyles, {'-','-.','-.','-.','--','--',':',':',':','-.','-.','--','--',':',':','-','-'}];
Colors = [ 204 0 0;
255 128 0;
0 153 153;
0 0 153;
127 0 255;
%51 255 51;
0 0 0;
153 153 0]/255;

Colors = [Colors ; Colors];
Colors = [Colors ; Colors];
Colors = [Colors ; Colors];

for i=1:min(numel(LineStyles),sum(relexch))
p(i).LineStyle = LineStyles{i};
p(i).Color = Colors(i,:);
end
caxes = gca;
xlabel('Oxygen uptake rate in ($\mu mol~g^{-1}~h^{-1}$)','Interpreter','latex','FontSize',20);
ylabel('Flux through transporter ($\mu mol~g^{-1}~h^{-1}$)','Interpreter','latex','FontSize',20);
hold on
set(caxes,'FontSize',18);
set(caxes,'Position',[0.1 0.19 0.87 0.74],'LineWidth',2,'TickLabelInterpreter','latex')
set(caxes,'XLim',[Vals(1) Vals(end)]);

leg = legend([RxnNames(relreacs)]);
leg.Interpreter = 'latex';
leg.Position = [0.67 0.709 0.3 0.2];


%% add hatched area and maximal literature value



hold on
set(caxes,'XLim',[floor(Vals(1)) maxscanval]);
optOxLine = plot([maxOxUptake maxOxUptake], caxes.YLim, 'LineWidth',3,'LineStyle',':','Color',[0 0 0]);
leg = legend([RxnNames(relreacs); 'Maximal Oxygen from Literature']);
leg.Interpreter = 'latex';
leg.Position = [0.61 0.709 0.25 0.2];
%s = surface([floor(Vals(1)) Vals(1)],caxes.YLim);

h1 = patch([floor(Vals(1)) Vals(1) Vals(1) floor(Vals(1))],[min(caxes.YLim), min(caxes.YLim), max(caxes.YLim), max(caxes.YLim)],'black');
hatch(h1, 45, [1 1 1], '-',3,1); 
%(OBJ,ANGLE,COLOR,STYLE,STEP,WIDTH)
%a1 =area(caxes, [floor(Vals(1)) Vals(1)], [max(caxes.YLim), max(caxes.YLim)], 'FaceColor','k','EdgeColor','none')
%a2 =area(caxes, [floor(Vals(1)) Vals(1)], [min(caxes.YLim), min(caxes.YLim)], 'FaceColor','k','EdgeColor','none')



%% Without Ala Dehydrogenase
SymModelWithoutAlaDehydrog = SymbioticModel;
AlaDH = ismember(SymModelWithoutAlaDehydrog.rxns,'ALANINE-DEHYDROGENASE-RXN_S');
SymbioticModel.ub(ismember(SymbioticModel.rxns,'TSyE_Poly-Hydroxybutyrate')) = 0;

SymModelWithoutAlaDehydrog.lb(AlaDH) = 0;
SymModelWithoutAlaDehydrog.ub(AlaDH) = 0;
SymModelWithoutAlaDehydrog.c(:) = 0;
SymModelWithoutAlaDehydrog.c(OxygenPos) = -1;
symsolWOAla = parsimOpt(SymModelWithoutAlaDehydrog);
SymModelWithoutAlaDehydrog.c(:) = 1;
symscansWOAla = ScanSingleReactionQP(SymModelWithoutAlaDehydrog,'TRSy_OXYGEN-MOLECULE',symsolWOAla.x(OxygenPos),maxscanval,scanpoints);


%% 
ReacIDs = {'TSyR_AMMONIA';'TSyR_L-ALPHA-ALANINE';'TRSy_VAL';'TRSy_LEU';'TRSy_ILE';'TRSy_OXYGEN-MOLECULE';'TRSy_GLT';'TSyR_L-ASPARTATE';'TSyR_AMMONIUM';'TRSy_MAL';'TRSy_SUC';'TSyE_Poly-Hydroxybutyrate'};
RxnNames = {'Ammonia export';'Alanine export';'Valine uptake';'Leucine uptake';'Isoleucine uptake';'Oxygen uptake';'Glutamate uptake';'Aspartate exchange';'Ammonium export';'Malate uptake';'Succinate uptake';'PHB production'};
AminoAcidUsage =  [findRxnsFromMets(SymbioticModel,{'GLT[S]','2-KETOGLUTARATE[S]'},'containsAll',1) ; findRxnsFromMets(SymbioticModel,{'L-ALPHA-ALANINE[S]'},'containsAll',1) ;  findRxnsFromMets(SymbioticModel,{'VAL[S]'},'containsAll',1)]; 

CarboxylateUsage = [findRxnsFromMets(SymbioticModel,{'OXALACETIC_ACID[S]'},'containsAll',1)];
%ReacIDs = [ReacIDs ; CarboxylateUsage];

%RxnNames = [RxnNames ; CarboxylateUsage];

exchangers  = ismember(SymbioticModel.rxns, ReacIDs);
%use the stds to see if there is any change in the values)
stds = mean(symscansWOAla')';
relexch = exchangers & abs(stds) > 1e-4;
relexch(OxygenPos) = 0;
[rpres,rpos] = ismember(SymbioticModel.rxns(relexch),ReacIDs);
relreacs = rpos(rpres);
%% Create Figure
%close all
f = figure;
f.Position = [50 50 1200 600];
Vals = linspace(symsolWOAla.x(OxygenPos),maxscanval,scanpoints);
%normscans = normalizeScans(symscanns);

Colors = [ 204 0 0;
255 128 0;
0 153 153;
0 0 153;
127 0 255;
%51 255 51;
0 0 0;
153 153 0]/255;

Colors = [Colors ; Colors];
Colors = [Colors ; Colors];
Colors = [Colors ; Colors];
p = plot(Vals,symscansWOAla(relexch,:),'LineWidth',3); 
LineStyles = {'-','-','-.','-.','--','--',':',':',':','-.','-.','--','--',':',':','-','-'};
for i=1:min(numel(LineStyles),sum(relexch))
p(i).LineStyle = LineStyles{i};
p(i).Color = Colors(i,:);
end
caxes = gca;
xlabel('Oxygen uptake rate in ($\mu mol~g^{-1}~h^{-1}$)','Interpreter','latex','FontSize',20);
ylabel('Flux through transporter ($\mu mol~g^{-1}~h^{-1}$)','Interpreter','latex','FontSize',20);
hold on
set(caxes,'FontSize',18);
set(caxes,'Position',[0.1 0.19 0.87 0.74],'LineWidth',2,'TickLabelInterpreter','latex')
set(caxes,'XLim',[Vals(1) Vals(end)]);

leg = legend([RxnNames(relreacs)]);
leg.Interpreter = 'latex';
leg.Position = [0.67 0.709 0.3 0.2];


%% add hatched area and maximal literature value



hold on
set(caxes,'XLim',[17 maxscanval]);
set(caxes,'YLim',[0 25]);
optOxLine = plot([maxOxUptake maxOxUptake], caxes.YLim, 'LineWidth',3,'LineStyle',':','Color',[0 0 0]);
leg = legend([RxnNames(relreacs); 'Maximal Oxygen from Literature']);
leg.Interpreter = 'latex';
leg.Position = [0.6 0.709 0.25 0.2];
%s = surface([floor(Vals(1)) Vals(1)],caxes.YLim);

h1 = patch([17 Vals(1) Vals(1) 17],[min(caxes.YLim), min(caxes.YLim), max(caxes.YLim), max(caxes.YLim)],'black');
hatch(h1, 45, [1 1 1], '-',3,1); 
%(OBJ,ANGLE,COLOR,STYLE,STEP,WIDTH)
%a1 =area(caxes, [floor(Vals(1)) Vals(1)], [max(caxes.YLim), max(caxes.YLim)], 'FaceColor','k','EdgeColor','none')
%a2 =area(caxes, [floor(Vals(1)) Vals(1)], [min(caxes.YLim), min(caxes.YLim)], 'FaceColor','k','EdgeColor','none')
