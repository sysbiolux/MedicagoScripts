function initSystem()

fastcorepathes = which('fastcore','-all')

for i=1:numel(fastcorepathes)    
if isempty(regexp(fastcorepathes{i},'/modelGeneration/fluxConsistency/FASTCORE'))
    usedpath = strrep(fastcorepathes{i},'fastcore.m','');
end
end
%Set the Fastcore path to the non Cobra version.
addpath(usedpath)