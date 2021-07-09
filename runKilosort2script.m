
%% add paths
marmoPath = addMarmoPipe();

user = 'gravedigger';
SERVER_DATA_DIR = addKiloPaths(user);
%% Select folder and convert raw data to a single binary file
if ~exist('DataFolder', 'var')
    assert(exist('SERVER_DATA_DIR', 'var')==1, 'Kilosort2ImportScript: SERVER_DATA_DIR does not exist. you must point to the raw data before running this')
    DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
end

%% if you need to test the channel map, use this cell


% Load up 10,000 samples of the channels
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
shanks2test = {'Nandy64', ...
    'Nandy', ...
    'Nandyflip', ...
    'NandyFrontFlip', ...
    'NeuroNexus64_2_32_flipped', ...
    'Nandy64headstage'};

for i = 1:numel(shanks2test)
   
    shank = hardware.electrodeFactory(shanks2test{i});
   
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


%% Once you know the channel map, import the data to single .dat file

% set the electrode to the one that worked in the cell above
shank = hardware.electrodeFactory('NandyFrontFlip'); 

io.oe2dat(DataFolder, shank, 'overwrite', true);
ops = io.loadOps(DataFolder);
SavePath = ops.root; % location of the processed data
currDir = pwd;
%% plot again the raw just to check that the channel map is correct
inds = (1:30e3) + 55000e3; % sample index
%%
% inds = inds + 10e3;
ops = ops(1);
inVolts = true;
data = io.loadRaw(ops, inds, inVolts, false);


% data = diff(data);
figure(1); clf
channelOffsets =(1:ops.Nchan)*2000;
channelOffsets = channelOffsets(1:size(data));
plot(bsxfun(@plus, data', channelOffsets))
%% Spike-sort single session

% DIRECTORY ON LOCAL DISK
fpath = DataFolder; %'C:\Raw\Logan_2020-03-04_09-51-54_neuronexus_D8\';
ops = io.loadOps(fpath);

% handle copying from path to path
ops = io.convertOpsToNewDirectory(ops, fpath);

%%  check the raw traces: look for dead channels, anything weird
[raw, timestamps] = io.loadRaw(ops, [1 10e3], true, false);

figure(1); clf
imagesc(raw)
title('Raw: look for dead channels')
xlabel('samples')
ylabel('channel')
numChan = size(raw,1);
v = var(raw, [],2);
%% find good channels and dead channels

figure(2); clf

load(ops.chanMap)
xcs = unique(xcoords);
for i = 1:numel(xcs)
    ix= xcoords == xcs(i);
    plot(mean(diff(raw(ix,:)),2)); hold on
    raw(ix,:) = raw(ix,:) - [0; mean(diff(raw(ix,:)),2)];
end
deadChannels = find(diff(v)>2000);
goodChannels = setdiff(1:numChan, deadChannels);
xlabel('Channel')

figure(1); clf
imagesc(raw); hold on
for i = 1:numel(deadChannels)
    plot(xlim, deadChannels(i)*[1 1], 'r')
end


figure(3); clf
plot(raw(deadChannels,:)')

%% plot again
clf
imagesc(raw(goodChannels,:)-median(raw(goodChannels,:)))

% plot(abs(zscore(var(raw, [],2)))>2)
%% modify ops struct for preprocessing
ops.deadChannels = deadChannels;
if isfield(ops, 'deadChannels')
    goodChannels = setdiff(1:ops.Nchan, ops.deadChannels);
end

channelOffsets =(1:ops.Nchan)*900;

% index
t = 2e3;
inds = (1:100e3) + 10000e3; % sample index
inds = inds + t;
inVolts = false;

% load raw
data = double(io.loadRaw(ops, inds, inVolts, false))';

ops.fslow = 7.5e3;
% ops = rmfield(ops, 'fslow');
ops.fshigh = 300;
ops.gfilter = 0; %[5 2]; %[5 1]; %[.1 1]; %[1 .1]; %[1 1]; %[.1 1];
datr = data(:,goodChannels);


% IF FILTER WITH JRCLUST
filtOpts.freqLimBP = [ops.fshigh ops.fslow];
filtOpts.freqLimNotch = [];
filtOpts.freqLimStop = [];
filtOpts.nSamplesPad = 100;
filtOpts.useGPUFilt = true;
filtOpts.sampleRate = ops.fs;
filtOpts.gainBoost = 1;
filtOpts.useElliptic = true;
filtOpts.filtOrder = 3;
% datr = filtfiltChain(datr, filtOpts);


% OR TRY DIFFERENTIATING FILTER FOM MATLAB
Fpass = 500;
Fstop = 1e3;
Nf = 8;
d = designfilt('differentiatorfir','FilterOrder',Nf, ...
    'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
    'SampleRate',ops.fs);

datr = filter(d, datr)*5;

% datr = gpufilter(data', ops, find(goodChannels));


datr = doFFTClean(datr, 5, false);

% DO LCAR AVERAGING
nchannel = 32;
nskip = 3;
nav = 8;
I = toeplitz([zeros(1,nskip) ones(1,nav) zeros(1,nchannel-nav-nskip)]);
I = blkdiag(I,I);
I = I(goodChannels,goodChannels);
I = I ./ sum(I);
datr = bsxfun(@minus, datr, datr*I); % common average reference

datr = zscore(datr)*50;

% norder = 2;
% datr = jrclust.filters.sgFilter(datr, norder)/(norder^2);

figure(1); clf

plot(datr + channelOffsets(goodChannels), 'k'); hold on

% datr = jrclust.utils.tryGpuArray(datr, true);
        
siteRMS = jrclust.utils.estimateRMS(datr, 10e5);
qqFactor = 4;
siteThresh = siteRMS*qqFactor;

inds = find(datr < -siteThresh);
[i, j] = ind2sub(size(datr), inds);
choff = channelOffsets(goodChannels);
plot(i, choff(j), 'r.')

nt = size(datr, 1);
nospks = setdiff(1:nt, unique(i + (-30:30)));

CC = cov(datr(nospks,:));
C0 = cov(datr);
figure(2); clf
subplot(1,3,1)
imagesc(CC)
title('spikes removed')
subplot(1,3,2)
imagesc(C0)
subplot(1,3,3)
imagesc(CC-C0)

chmap = load(ops.chanMap);
yc = chmap.ycoords(goodChannels);
xc = chmap.xcoords(goodChannels);
nRange = 14;
Wrot = whiteningLocal(gather(CC), yc, xc, nRange);
for i = 1:size(Wrot,2)
    Wrot(:,i) = Wrot(:,i)/norm(Wrot(:,i));
end
% Wrot = Wrot*50;

figure(3); clf
subplot(2,1,1)
imagesc((datr)');
xlim([0 500])
subplot(2,1,2)

datr = datr*Wrot;
imagesc((datr)')
xlim([0 500])
colormap jet

figure(1)
plot(datr + channelOffsets(goodChannels), 'b'); hold on


%%
% ops.deadChannels = deadChannels;
pathToYourConfigFile = 'C:\Users\Jake\Documents\MATLAB\MarmoPipe\KiloSort2\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'configMitchellLab.m'))
ops.NchanTOT = ops.Nchan;
ops.trange = [1 inf];

%%
% % If you already ran Kilosort, you can load rez and run from here
% if exist(fullfile(ops.root, 'rez.mat'), 'file') && exist(ops.fproc, 'file')
%     load(fullfile(ops.root, 'rez.mat'))
% else
%     % preprocess data to create temp_wh.dat
%     rez = preprocessDataSub(ops);
% 
%     % time-reordering as a function of drift
%     rez = clusterSingleBatches(rez);
%     save(fullfile(ops.root, 'rez.mat'), 'rez', '-v7.3');
%     
% end

%%
ops.trange      = [0 Inf]; % TIME RANGE IN SECONDS TO PROCESS
 % these settings overwrite any settings from config
ops.lam      = 10;   % weighting on the amplitude penalty (like in Kilosort1, but it has to be much larger)
ops.Th      = [10 4];  % lower bound on acceptable single spike quality
% ops.ThS      = [6 6];  % lower bound on acceptable single spike quality
ops.deadChannels = [];
ops.momentum = [50 400]; % number of samples to average over
ops.Nfilt = 2*ops.Nchan;
ops.minFR = 1;
ops.minfr_goodchannels = 0.01;
ops.nSkipCov = 100;
ops.fslow = 5500;
% ops = rmfield(ops, 'fslow');
ops.scaleproc = 50;
ops.fshigh = 300;
ops.artifactThresh = 250;
ops.artifactNchans = 20;
ops.debug = false;
ops.whiteningRange = 8; % negative means custom rotation matrix
ops.reorder = 1;
ops.gfilter = 0;

%% preprocess
rez = preprocessDataSub(ops);

%%
% time-reordering as a function of drift
rez = clusterSingleBatches(rez);
save(fullfile(ops.root, 'rez.mat'), 'rez', '-v7.3');
ibatch = 1;
disp('Done clustering and saving rez')
%% look at the filtered/whitened data
ibatch = ibatch + 1; % step through batches
dat = get_batch(rez.ops, ibatch);
figure(1); clf
imagesc(dat')
xlim([0 2e3])
colormap jet % use a nonlinear colormap to see what thresholding will do
% plot(dat)
%%  main sorting
load(fullfile(ops.root, 'rez.mat'))
rez.ops.Th      = [12 4];  % lower bound on acceptable single spike quality
rez.ops.lam = 2;
ops.momentum = [20 400];
% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

%% run matching pursuit 

% rez.ops.Th      = [12 3];  % lower bound on acceptable single spike quality
% rez.ops.lam = 10;
% rez.ops.momentum = [20 400];
% rez = runTemplates(rez);


%% plot waveforms templates
Nch = numel(rez.xc);
choffsets = (1:Nch)*5;
cids = unique(rez.st3(:,2));
NC = numel(cids);
cmap = lines(NC);
figure(1); clf
for cc = 1:NC
    wf = rez.Wraw(:,:,cids(cc));
    xax = linspace(-1/2, 1/2, size(rez.Wraw, 2)) + cc;
    plot(xax, wf' + choffsets, 'Color', cmap(cc,:)); hold on
    text(cc, choffsets(end) + 10, sprintf('%d', cc), 'Color', cmap(cc,:))
end
%% plot autocorrelation and ccgs
figure(2); clf
cids = unique(rez.st3(:,2));
NC = numel(cids);
nbins = 100;
fs = 30e3;
cmap = lines(NC);
for ii = 1:NC
    fprintf('%d/%d\n', ii, NC)
    sp1 = rez.st3(rez.st3(:,2) == cids(ii))/fs;
    if numel(sp1) ==0
        continue
    end
    
    for jj = 1:NC
        
        if jj < ii
            continue
        end
        
        sp2 = rez.st3(rez.st3(:,2) == cids(jj))/fs;
        
        if numel(sp2) ==0
            continue
        end
        K = ccg(sp1, sp2, nbins, 1e-3);
        if ii == jj
            K(nbins + 1) = 0;
        end
        K = (K - min(K)) ./ (max(K) - min(K));
        xax = linspace(-1/2, 1/2, nbins*2 + 1) + jj;
        
        plot(xax, K-ii, 'Linewidth', 1, 'Color', cmap(jj,:)); hold on%'BarWidth', 1, 'EdgeColor', 'None', 'FaceColor', cmap(jj,:)); hold on
        
        
    end
end
%% final merges
rez0 = rez;
%%
rez = find_merges_jly(rez, 1, true);

%% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

%% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))
fprintf('found %d total units\n', numel(unique(rez.st3(:,2))))
% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, ops.root);

%%
fclose all
delete(ops.fproc)
ops.fshigh = 300;
preprocess.save_filtered_data(ops)
fprintf('Done\n')
%%
rezToPhy(rez, ops.root);
%%
fprintf('Once you finish sorting in the GUI and you have saved the output\n')
fprintf('Open Anaconda Prompt and enter\n')
fprintf('activate phy2\n')
fprintf('cd %s\n', ops.root)
fprintf('phy template-gui params.py\n')
% 
%%
sp = io.getSpikesFromKilo(ops);

%% plot raw data
[ss, inds] = sort(rez.st3(:,1));
clu = rez.st3(inds,2);
cids = unique(clu);
NC = numel(cids);

batchinds = [1 10e3];
data = io.loadRaw(ops, batchinds, true, false);

data = data';
figure(10); clf
cmap = lines;
cmap = zeros(size(cmap));
for i = 5
[b,a] = butter(i, 300/30000*2, 'high');
data = filter(b, a, data);

plot(bsxfun(@plus, data, (1:size(data,2))*200), 'Color', cmap(i,:)); hold on
end

%%

ix = ss > batchinds(1) & ss < batchinds(2);
t0 = batchinds(1)-1;

win = -10:30;
nsamp = numel(win);
nch = size(data,2);
wfs = zeros(nsamp, nch, NC);

cmap = lines;
for cc = 1:NC
    spix = ss(clu==cids(cc));
    plot((spix*[1 1])', ylim, 'Color', cmap(cc,:))
    spix(spix<11) = [];
    spix(spix> (batchinds(2) -30)) = [];
    nspk = numel(spix);
    for ch = 1:nch
        wfs(:,ch,cc) =  mean(reshape(data(spix(:) + win,ch), [nspk, nsamp]));
    end
    drawnow
end
    

figure(3); clf
for cc= 1:NC
    plot(wfs(:,:,cc));
    pause
end
% plot(ss(ix)-t0


plot(ss, rez.st3(inds,2), '.')

