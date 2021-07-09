ops = io.loadOps(DataFolder);
% overwrite ops
pathToYourConfigFile = 'ExternalRepos\KiloSort2\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'configMitchellLab.m'))
ops.NchanTOT = ops.Nchan;
ops.trange = [1 inf];

% If you already ran Kilosort, you can load rez and run from here
if exist(fullfile(ops.root, 'rez.mat'), 'file') && exist(ops.fproc, 'file')
    load(fullfile(ops.root, 'rez.mat'))
else
    % preprocess data to create temp_wh.dat
    rez = preprocessDataSub(ops);

    % time-reordering as a function of drift
    rez = clusterSingleBatches(rez);
    save(fullfile(ops.root, 'rez.mat'), 'rez', '-v7.3');
    
end

%%
% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

%% final merges
rez = find_merges(rez, 1);

%% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, ops.root);

%% if you want to save the results to a Matlab file... 

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(ops.root, 'rez2.mat');
save(fname, 'rez', '-v7.3');

