
marmoPath = addMarmoPipe();

user = 'gravedigger';
SERVER_DATA_DIR = addKiloPaths(user);
%% Select folder and convert raw data to a single binary file
if ~exist('DataFolder', 'var')
    assert(exist('SERVER_DATA_DIR', 'var')==1, 'Kilosort2ImportScript: SERVER_DATA_DIR does not exist. you must point to the raw data before running this')
    DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
end

shank = hardware.electrodeFactory('Nandy64headstage'); 

io.oe2dat(DataFolder, shank, 'overwrite', false);
ops = io.loadOps(DataFolder);
SavePath = ops.root; % location of the processed data
currDir = pwd;
%% plot the raw just to check that the channel map is correct
inds = (1:30e3) + 55000e3; % sample index
% inds = inds + 10e3;
ops = ops(1);
inVolts = true;
data = io.loadRaw(ops, inds, inVolts);
data = bsxfun(@minus, data, mean(data)); % common average reference
figure(1); clf
channelOffsets =(1:ops.Nchan)*200;
plot(bsxfun(@plus, data', channelOffsets))
% title(SavePath)

%% SPIKE SORTING - Do either a or b - DONT DO BOTH
%% a)  run Kilosort without GUI
mitchell_kilosort

%% b) run Kilosort

kiloDir = fileparts(which('kilosort'));
% needs to run from Kilosort directory?
cd(kiloDir)

fprintf('Opening Kilosort GUI\n')

kilosort

%% manual curation in phy
% 
fprintf('Once you finish sorting in the GUI and you have saved the output\n')
fprintf('Open Anaconda Prompt and enter\n')
fprintf('activate phy2\n')
fprintf('cd %s\n', SavePath)
fprintf('phy template-gui params.py\n')
% 
%% save spikes as struct
ops = io.loadOps(DataFolder);
info = io.loadEphysInfo(DataFolder);
sp = io.getSpikesFromKilo(ops, info);
% sp = loadKSdir(SavePath);
save(fullfile(SavePath, 'spikes-kilo.mat'), '-v7.3', '-struct', 'sp')

% 
% cd(currDir)
% sp = loadKSdir(SavePath);
% save(fullfile(SavePath, 'spikes-kilo.mat'), '-v7.3', '-struct', 'sp')
