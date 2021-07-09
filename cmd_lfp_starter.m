
marmoPath = addMarmoPipe();

user = 'gravedigger';
SERVER_DATA_DIR = addKiloPaths(user);
%% Select folder and convert raw data to a single binary file
if ~exist('DataFolder', 'var')
    assert(exist('SERVER_DATA_DIR', 'var')==1, 'Kilosort2ImportScript: SERVER_DATA_DIR does not exist. you must point to the raw data before running this')
    DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
end

ops = io.loadOps(DataFolder);
SavePath = ops.root; % location of the processed data
currDir = pwd;

%% get EXP struct

user = 'jakegravedigger';

switch user
    case 'judehome'
        SERVER_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Raw';
        PROCESSED_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Processed';
    case 'shannawork'
        SERVER_DATA_DIR = 'C:\Users\Shanna\Dropbox\Marmo Lab Website\PSA\DDPI_Raw';
        PROCESSED_DATA_DIR = 'C:\Users\Shanna\Dropbox\Marmo Lab Website\PSA\DDPI_Processed';
    case 'jakegravedigger'
        SERVER_DATA_DIR = 'C:\Raw';
        PROCESSED_DATA_DIR = 'Z:\Data\Processed'; %'C:\Processed';
end

% DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
[~, FileTag] = fileparts(DataFolder);
 DataFolder = fullfile(SERVER_DATA_DIR,FileTag);

ExpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'.mat'];

load(ExpFile);
%% plot the raw just to check that the channel map is correct
inds = (1:30e3) + 55000e3; % sample index

inVolts = true;
data = io.loadRaw(ops, inds, inVolts);

%% Get LFP
%[lfp, ~, lfpinfo] = io.getLFP(ops, false, false);
[lfp, ~, lfpinfo] = io.getLFP(ops, true, false);


%% Run 
CSD_Spike_Onset_Script.m 
%% Compute CSD
events = sort(reshape(CSD.Onsets_Ephys, 1, size(CSD.Onsets_Ephys,1)*size(CSD.Onsets_Ephys,2)));
window = [-100 200]; % window for computing CSD
%stats = csd.csdFlash.csdBasic(lfp(:,[33:33+10 33+12:64]), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);
deadCh = [3 7 11 15 19];
a = 1:32;
a(deadCh) = [];
%stats = csd.csdFlash.csdBasic(lfp(:,a), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);

%stats2 = csd.csdFlash.csdBasic(lfp(:,[33:39 41:43 45:64]), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);

stats = csd.csdFlash.csdBasic(lfp(:,1:32), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);
stats2 = csd.csdFlash.csdBasic(lfp(:,33:64), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);


%%

events = sort(reshape(CSD.Onsets_Ephys, 1, size(CSD.Onsets_Ephys,1)*size(CSD.Onsets_Ephys,2)));
window = [-100 200]; % window for computing CSD
%stats = csd.csdFlash.csdBasic(lfp(:,[33:33+10 33+12:64]), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);
deadCh = [3 7 11 15 19];
a = 1:32;
a(deadCh) = [];
%stats = csd.csdFlash.csdBasic(lfp(:,a), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);


fragments = round(lfpinfo.fragments);
ev = io.convertTimeToSamples(events, lfpinfo.sampleRate, lfpinfo.timestamps(:), fragments(:));

ix = nan(1,length(events));
for i = 1:length(events)
    [ ~, ix(i) ] = min( abs( lfp.timestamps-events(i) ) );
end


[~,~,lipo] = histcounts(events, lfp.timestamps);


%stats2 = csd.csdFlash.csdBasic(lfp(:,[33:39 41:43 45:64]), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);

%stats = csd.csdFlash.csdBasic(data1(:,1:32), events, lfpinfo, 'method' , 'spline', 'window', window, 'plot', false);

%% Plot CSD
figure(3)
subplot(1,2,1)
imagesc(stats.time, stats.chDepths, stats.CSD-mean(stats.CSD(:))); axis xy
colormap jet
hold on
plot(stats.time, bsxfun(@plus, stats.STA', stats.chDepths(1:size(stats.STA,1))), 'Color', repmat(.1, 1, 3))
xlim(window)
% plot(stats.time([1 end]), stats.sinkDepth*[1 1], 'w--', 'Linewidth', 1)
% plot(stats.time([1 end]), [1; 1]*stats.reversalPointDepth, 'r--', 'Linewidth', 1)
tmp = abs(stats.reversalPointDepth - stats.sinkDepth);
tmp = tmp + stats.sinkDepth;
colorbar

subplot(1,2,2)
imagesc(stats2.time, stats2.chDepths, stats2.CSD-mean(stats2.CSD(:))); axis xy
colormap jet
hold on
plot(stats2.time, bsxfun(@plus, stats2.STA', stats2.chDepths(1:size(stats2.STA,1))), 'Color', repmat(.1, 1, 3))
xlim(window)
% plot(stats.time([1 end]), stats.sinkDepth*[1 1], 'w--', 'Linewidth', 1)
% plot(stats.time([1 end]), [1; 1]*stats.reversalPointDepth, 'r--', 'Linewidth', 1)
tmp = abs(stats2.reversalPointDepth - stats2.sinkDepth);
tmp = tmp + stats2.sinkDepth;
colorbar
% plot(stats.time([1 end]), [1; 1]*tmp, 'r--', 'Linewidth', 1);

%% Compute band power

p = bandpower(stats.STA',info.sampleRate,[20 60]);
figure(2)
plot(1:size(stats.STA, 1),p)



