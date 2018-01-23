function MedicagoNew = importMedicago()
%Read the Medicago Model and bring it into a nicer Form.

MedicFromSBML = readCbModel(['Data' filesep 'MedicagoTruncatula.xml']);

Medic = MedicFromSBML;

%We have a xml useable format where any non sBML characters are replaced by
%their respective integer value.
origstrings = {'__45__','__46__','__40__','__41__','__43__'};
targetstrings = {'-','.','(',')','+'};
for i=1:numel(targetstrings)
    Medic.rxns = regexprep(Medic.rxns,origstrings{i},targetstrings{i});
    Medic.mets= regexprep(Medic.mets,origstrings{i},targetstrings{i});
end

Medic.rxns = regexprep(Medic.rxns,'^R_(T[A-Z][A-Z])','$1'); %replace all R_ before Transport Reactions.

Medic.rxns = regexprep(Medic.rxns,'^_','');
Medic.mets = regexprep(Medic.mets,'^_','');
%Fix the bounds (new Cobra version does not read kinetic parameters if FBC
%Version 2
Medic.lb(Medic.lb == -1000) = -99999;
Medic.ub(Medic.ub == 1000) = 99999;
Medic.osense = -1;

%Replace the name of the ATPase, which is originally RXN-11135_C
atpasepos = find(ismember(Medic.rxns,'RXN-11135_C'));
Medic.rxns{atpasepos} = 'ATPase';
Medic.lb(atpasepos) = 0;
%Also fix the DEHOG to zero flux (again from Kinetic parameters not read
%any more.
dehogpos = find(ismember(Medic.rxns,'DEHOG'));
Medic.lb(dehogpos) = 0;
Medic.ub(dehogpos) = 0;


%And change the BiomassWithOutStarch to BiomassShootWithOutStarch
biomasswostarchshootpos = find(ismember(Medic.rxns,'BiomassWithOutStarch'));
Medic.rxns{biomasswostarchshootpos} = 'BiomassShootWithOutStarch';
Biomassreacs = {'BiomassShoot','BiomassRoot', 'BiomassShootWithOutStarch'};
Medic.ub(ismember(Medic.rxns,Biomassreacs)) = 1000;
%Reset the objective...
Medic.c(:) = 0;
Medic.c(1) = 1; % The Ammonium importer.
MedicagoNew = Medic;
