
%%

marmoPath = addMarmoPipe();

user = 'gravedigger';
SERVER_DATA_DIR = addKiloPaths(user);
%% Select folder and convert raw data to a single binary file
if ~exist('DataFolder', 'var')
    assert(exist('SERVER_DATA_DIR', 'var')==1, 'Kilosort2ImportScript: SERVER_DATA_DIR does not exist. you must point to the raw data before running this')
    DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
end


%% Load up 10,000 samples of the channels
Nchan = 64;
nSamples = 20e3;
data = zeros(Nchan, nSamples);

for i = 1:Nchan
    fprintf('Loading channel %d\n', i)    
    chFiles = dir(fullfile(DataFolder, sprintf('*CH%d.continuous', i)));
    if ~isempty(chFiles)
        filename = fullfile(DataFolder, chFiles.name);
        
        data_seg = load_open_ephys_data_faster(filename, 'unscaledInt16');
        data(i,:) = double(data_seg(1:nSamples));
    end
end
    
%% Test channel map
shanks2test = {'Nandy', ...
    'Nandyflip', ...
    'NandyFrontFlip', ...
    'NeuroNexus64_2_32_flipped', ...
    'Nandy64headstage'};

shanks2test = {'NeuroNexus64_2_32', ...
        'Nandy', ...
        'Nandy64', ...
        'Nandyflip', ...
        'NandyFrontFlip', ...
        'NeuroNexus64_2_32_flipped' ...
        };

for i = 1:numel(shanks2test)
    % shank = hardware.electrodeFactory('Nandy64');
    shank = hardware.electrodeFactory(shanks2test{i});
    % shank = hardware.electrodeFactory('Nandyflip');
%     shank{1}.headstages = {hardware.headstage.intan_RHD2164FLIP};
    headstage = hardware.combineHeadstages(shank{1}.headstages);
    
    channelMap = headstage.mapChannels(shank{1});
    
    sprintf('%d,', channelMap)
    
    figure(i); clf
    channelOffsets =(1:Nchan)*1000;
    dtmp = data(channelMap,:)';
    dtmp = dtmp - mean(dtmp,2);
    plot(bsxfun(@plus, dtmp, channelOffsets))
    title(strrep(shanks2test{i}, '_', '-'))
end

% plot(bsxfun(@plus, double(squeeze(DATA(:,:,1))), 10*channelOffsets(1:61)))
% [map.chanMap(:) map.connected(:)]

%% Test channel map - GABE
%     'Nandy', ...
%         'Nandyflip', ...
%         'NandyFrontFlip', ...
%         'NeuroNexus64_2_32_flipped'
% shank = hardware.electrodeFactory('NandyFrontFlip');
shank = cell(1);
shank{1} = hardware.electrode.NeuroNexus64_2_32_frontflip;
shank{1}.headstages{1} = hardware.headstage.intan_RHD2132;

headstage = hardware.combineHeadstages(shank{1}.headstages);

channelMap = headstage.mapChannels(shank{1});

sprintf('%d,', channelMap)

figure(1); clf
subplot(1,2,1)
channelOffsets =(1:Nchan)*2000;
dtmp = data(channelMap,:)';
dtmp = dtmp - mean(dtmp,2);
plot(bsxfun(@plus, dtmp, channelOffsets))
title('NandyFrontFlip')


shank = hardware.electrodeFactory('Nandyflip');

headstage = hardware.combineHeadstages(shank{1}.headstages);

channelMap = headstage.mapChannels(shank{1});

sprintf('%d,', channelMap)

figure(1);
subplot(1,2,2)
channelOffsets =(1:Nchan)*2000;
dtmp = data(channelMap,:)';
dtmp = dtmp - mean(dtmp,2);
plot(bsxfun(@plus, dtmp, channelOffsets))
title('Nandyflip')



% plot(bsxfun(@plus, double(squeeze(DATA(:,:,1))), 10*channelOffsets(1:61)))
% [map.chanMap(:) map.connected(:)]

%% save channel map to file
% offset channels > 32
channelMap(channelMap > 32) = channelMap(channelMap > 32) + 3;
%%


for i = 1:numel(channelMap)
    fprintf('%d) %d\n', i, channelMap(i))
end

