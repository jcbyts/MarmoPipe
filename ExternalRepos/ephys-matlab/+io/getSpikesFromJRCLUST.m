function sp = getSpikesFromJRCLUST(ops, info)
% getSpikesFromKilo gets the spike sorting output from Kilosort
% requires the spikeTools from cortexlab github

if isa(ops, 'table') %hold over from when Jake tracked everything with a .csv file
    assert(any(strcmp(ops.SpikeSorting, 'JRC')), 'JRCLUST has not been run on this session')
    ops = io.loadOps(ops);
    sp = [];
    for i = 1:numel(ops)
        try
            sp = [sp get_spikes_jrc_helper(ops(i))];
        end
    end
    
elseif isstruct(ops)
    if ~exist('info', 'var')
        info = load(fullfile(ops.root, 'ephys_info.mat'));
    end
    sp = get_spikes_jrc_helper(ops, info);
end


function sp = get_spikes_jrc_helper(ops, info)


flist = dir(fullfile(ops.root, '*_res.mat'));
fnames = unique(arrayfun(@(x) x.name(1:strfind(x.name, '_')-1), flist, 'uni', 0));

[~, fid] = max([flist.datenum]);
fname = fullfile(ops.root, sprintf('%s_res.mat', fnames{fid}));
assert(exist(fname, 'file')==2, 'results file does not exist')

rez = load(fname);

if nargin < 2
    info = load(fullfile(ops.root, 'ephys_info.mat'));
end


folderNames = {ops.root};

% get configuration info
cfg = jrclust.utils.mToStruct(fullfile(ops.root, [fnames{1} '.prm']));
Fs = cfg.sampleRate;

% channel geometry
xc = cfg.siteLoc(:,1);
yc = cfg.siteLoc(:,2);

for f = 1:numel(folderNames)    

cids = (1:numel(rez.clusterNotes))';
nClusters = numel(cids);
cgs = 3*ones(nClusters, 1); % 0 = noise; 1 = MUA; 2 = Good; 3 = Unsorted
labels = {'noise', 'multi', 'single'};
for i = 1:numel(labels)
    ix = strcmp(rez.clusterNotes, labels{i});
    cgs(ix) = i-1;
end

good = cgs > 0 & cgs < 3;
cgs = cgs(good);
cids = cids(good);
nClusters = numel(cids);

% ss is a length nSpikes vector with the spike time of every spike (in
% samples)
ss = [];
% clu is a length nSpikes vector with the cluster identities of every
% spike
clu = [];
for cc = 1:nClusters
   ss = [ss; double(rez.spikeTimes(rez.spikesByCluster{cids(cc)}))]; %#ok<*AGROW>
   clu = [clu; double(rez.spikeClusters(rez.spikesByCluster{cids(cc)}))];
end

[ss, ind] = sort(ss);
clu = clu(ind);

% convert to times in seconds
st = io.convertSamplesToTime(ss, Fs, info.timestamps(:), info.fragments(:));

% 0 = noise; 1 = MUA; 2 = Good; 3 = Unsorted

% Noise spikes:
% should be excluded from everything; we do this in a moment.

% Both MUA and Unsorted:
% reflect real spikes (in my judgment) but ones that
% couldn't be isolated to a particular neuron. They could include in
% analyses of population activity, but might include spikes from multiple
% neurons (or only partial spikes of a single neuron. 

% Good clusters:
% ones I judged to be well-isolated based on a combination of subjective
% criteria: how clean the refractory period appeared, how large the spike
% amplitudes were, how unique the waveform shapes were given the
% surrounding context.
    

savefname = 'spJRC.mat';
    
clusterDepths = rez.clusterCentroids(cids,2);
clusterX = rez.clusterCentroids(cids,1);

spikeAmps = cfg.bitScaling * double(rez.spikeAmps);
spikeDepths = rez.spikePositions(:,2)';
spikeX = rez.spikePositions(:,1)';

% remove noise clusters
sp(f).name = ops.root;
sp(f).sorter = 'JRCLUST';
sp(f).clu = clu;
sp(f).ss = ss;
sp(f).st = st;
sp(f).spikeTemplates = clu;
sp(f).cgs = cgs;
sp(f).cids = cids;
sp(f).yc = yc;
sp(f).xc = xc;
sp(f).ycoords = yc;
sp(f).xcoords = xc;
sp(f).spikeAmps = spikeAmps;
sp(f).spikeDepths = spikeDepths;
sp(f).spikeX = spikeX;
sp(f).tempsUnW = permute(rez.meanWfGlobal(:,:,cids), [3 1 2]);

sp(f).cgs2 = cgs; % duplicate the rating
sp(f).cids2 = cids;
sp(f).uQ = rez.unitSNR(cids);
sp(f).cR = rez.unitIsoDist(cids);
sp(f).isiV = rez.unitISIRatio(cids);

sp(f).firingRates = rez.unitCount(cids) / (max(st)-min(st));
sp(f).clusterDepths = clusterDepths;
sp(f).clusterX = clusterX;

end


%% save
for s = 1:numel(folderNames)
    sptmp = sp(s);
    save(fullfile(folderNames{s}, savefname), '-v7.3', '-struct', 'sptmp');
end

