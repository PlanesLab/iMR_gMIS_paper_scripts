
nLayers = 1;
name = 'Human-GEM-1.14.0'
numWorkers = 8;

modelName1 = 'TRRUST'
% modelName1 = 'Omnipath'
% modelName1 = 'Dorothea'
% modelName1 = 'Signor'
signalingNetwork = readtable(['./Databases/' modelName1 '.txt']);
load([name '.mat']);

modelName = [modelName1 '_' name '_' type]

type = 'SL';
% type = 'SDL'


maxlength = 5;

if strcmp(type, 'SDL') 
    [G, G_ind, G_rxns, ~, ~, G_time, ~] = buildGmatrix_iMRmodel_gMIS(modelName, model, signalingNetwork, nLayers, 774, numWorkers);
elseif strcmp(type, 'SL')
    [G, G_ind, G_rxns, ~, ~, G_time, ~] = buildGmatrix_iMRmodel_gMCS(modelName, model, signalingNetwork, nLayers, 774, numWorkers);
end

modelName_wlayer = [modelName '_' num2str(nLayers) '_layers']
filterGmatrix(G, G_ind, maxlength, G_time, G_rxns, modelName_wlayer);

[gmcs, gmcs_time, gmcs_onlynew] = calculateGeneMCS(modelName_wlayer, model, 5000, 10, 'timelimit', 300, 'numWorkers', numWorkers, 'forceLength', 0);

save(['Results' filesep modelName_wlayer '.mat'],'model', 'gmcs', 'gmcs_time', 'gmcs_onlynew');

