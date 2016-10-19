function Medicago = importMedicago()
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

Medic.rxns = regexprep(Medic.rxns,'^_','');
Medic.mets = regexprep(Medic.mets,'^_','');
%Replace the name of the ATPase, which is originally RXN-11135_C
atpasepos = find(ismember(Medic.rxns,'RXN-11135_C'));
Medic.rxns{atpasepos} = 'ATPase';

%And change the BiomassWithOutStarch to BiomassShootWithOutStarch
biomasswostarchshootpos = find(ismember(Medic.rxns,'BiomassWithOutStarch'));
Medic.rxns{biomasswostarchshootpos} = 'BiomassShootWithOutStarch';

Medicago = Medic;
