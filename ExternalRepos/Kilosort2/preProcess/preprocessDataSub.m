function rez = preprocessDataSub(ops)
% this function takes an ops struct, which contains all the Kilosort2 settings and file paths
% and creates a new binary file of preprocessed data, logging new variables into rez.
% The following steps are applied:
% 1) conversion to float32;
% 2) common median subtraction;
% 3) bandpass filtering;
% 4) channel whitening;
% 5) scaling to int16 values

tic;
ops.nt0 	  = getOr(ops, {'nt0'}, 61); % number of time samples for the templates (has to be <=81 due to GPU shared memory)
ops.nt0min  = getOr(ops, 'nt0min', ceil(20 * ops.nt0/61)); % time sample where the negative peak should be aligned

NT       = ops.NT ; % number of timepoints per batch
NchanTOT = ops.NchanTOT; % total number of channels in the raw binary file, including dead, auxiliary etc

bytes       = get_file_size(ops.fbinary); % size in bytes of raw binary
nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
ops.tstart  = ceil(ops.trange(1) * ops.fs); % starting timepoint for processing data segment
ops.tend    = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); % ending timepoint
ops.sampsToRead = ops.tend-ops.tstart; % total number of samples to read
ops.twind = ops.tstart * NchanTOT*2; % skip this many bytes at the start

Nbatch      = ceil(ops.sampsToRead /(NT-ops.ntbuff)); % number of data batches
ops.Nbatch = Nbatch;

[chanMap, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(ops.chanMap); % function to load channel map file
ops.NchanTOT = getOr(ops, 'NchanTOT', NchanTOTdefault); % if NchanTOT was left empty, then overwrite with the default

% first find dead channels
if ~isfield(ops, 'deadChannels')
    igood = get_dead_channels(ops, chanMap);
    chanMap = chanMap(igood);
else
    chanMap = setdiff(chanMap, ops.deadChannels);
end
% igood = get_dead_channels(ops, chanMap);
% chanMap = chanMap(igood);

if getOr(ops, 'minfr_goodchannels', .1)>0 % discard channels that have very few spikes
    
    % determine bad channels
    fprintf('Time %3.0fs. Determining good channels by firing rate.. \n', toc);
    igood = get_good_channels(ops, chanMap);

    chanMap = chanMap(igood); %it's enough to remove bad channels from the channel map, which treats them as if they are dead

    xc = xc(igood); % removes coordinates of bad channels
    yc = yc(igood);
    kcoords = kcoords(igood);

    ops.igood = igood;
else
    ops.igood = true(size(chanMap));
end

ops.Nchan = numel(chanMap); % total number of good channels that we will spike sort
ops.Nfilt = getOr(ops, 'nfilt_factor', 4) * ops.Nchan; % upper bound on the number of templates we can have

rez.ops         = ops; % memorize ops

rez.xc = xc; % for historical reasons, make all these copies of the channel coordinates
rez.yc = yc;
rez.xcoords = xc;
rez.ycoords = yc;
% rez.connected   = connected;
rez.ops.chanMap = chanMap;
rez.ops.kcoords = kcoords;


NTbuff      = NT + 4*ops.ntbuff; % we need buffers on both sides for filtering

rez.ops.Nbatch = Nbatch;
rez.ops.NTbuff = NTbuff;
rez.ops.chanMap = chanMap;


fprintf('Time %3.0fs. Computing whitening matrix.. \n', toc);

% this requires removing bad channels first

Wrot = get_whitening_matrix(rez); % outputs a rotation matrix (Nchan by Nchan) which whitens the zero-timelag covariance of the data

% override whitening with custom matrix
% tvec = [10 2 -1 zeros(1,numel(chanMap)-3)];
% tvec = tvec / norm(tvec);
% Wrot = toeplitz(tvec, tvec);

% jake custom: whitening rotates only -- normalize to unit vectors
% for i = 1:size(Wrot,2)
%     Wrot(:,i) = Wrot(:,i)/norm(Wrot(:,i));
% end

% JAKE CUSTOM FILTERING
Fpass = 500;
Fstop = 1e3;
Nf = 8;
dfilt = designfilt('differentiatorfir','FilterOrder',Nf, ...
    'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
    'SampleRate',ops.fs);

% % filtering from JRCLUST
% filtOpts.freqLimBP = [ops.fshigh ops.fslow];
% filtOpts.freqLimNotch = [];
% filtOpts.freqLimStop = [];
% filtOpts.nSamplesPad = 100;
% filtOpts.useGPUFilt = true;
% filtOpts.sampleRate = ops.fs;
% filtOpts.gainBoost = 1;
% filtOpts.useElliptic = true;
% filtOpts.filtOrder = 3;

% Local common average referencing
nchannel = 32;
nskip = 3;
nav = 8;
Icar = toeplitz([zeros(1,nskip) ones(1,nav) zeros(1,nchannel-nav-nskip)]);
if ops.Nchan > 32
    Icar = blkdiag(Icar,Icar);
elseif ops.Nchan==32
    % do nothing
else
    error('wrong number of channels. need to implement local common average referencing')
%     Icar = bldiag(Icar, Icar, Icar);
end

Icar = Icar(chanMap,chanMap);
Icar = Icar ./ sum(Icar); % take mean (excluding padding around focal channel)

fprintf('Time %3.0fs. Loading raw data and applying filters... \n', toc);

fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
fidW        = fopen(ops.fproc,   'w'); % open for writing processed data
debug = getOr(ops, 'debug', 0);

for ibatch = 1:Nbatch
    % we'll create a binary file of batches of NT samples, which overlap consecutively on ops.ntbuff samples
    % in addition to that, we'll read another ops.ntbuff samples from before and after, to have as buffers for filtering
    offset = max(0, ops.twind + 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff)); % number of samples to start reading at.
    if offset==0
        ioffset = 0; % The very first batch has no pre-buffer, and has to be treated separately
    else
        ioffset = ops.ntbuff;
    end
    fseek(fid, offset, 'bof'); % fseek to batch start in raw file

    buff = fread(fid, [NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
    if isempty(buff)
        break; % this shouldn't really happen, unless we counted data batches wrong
    end
    nsampcurr = size(buff,2); % how many time samples the current batch has
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr); % pad with zeros, if this is the last batch
    end

%     datr    = gpufilter(buff, ops, chanMap); % apply filters and median subtraction
    
%     datr = filtfiltChain(single(buff(chanMap,:))', filtOpts);

    datr = filter(dfilt, single(buff(chanMap,:))'); % CUSTOM FILTER
    
    datr    = datr(ioffset + (1:NT),:); % remove timepoints used as buffers
    
    
    datr = bsxfun(@minus, datr, datr*Icar); % common average reference
    
    
    
    if debug
        figure(1); clf
        d = datr;
        clim = [min(min(d)) max(max(d))];
        subplot(211)
        imagesc(d', clim)
        xlim([0 1e3])
%         plot(d + (1:size(d,2))*1300)
%         hold on
    end
    
    datr = doFFTClean(datr, 6, true);
%     norder = 2;
%     datr = jrclust.filters.sgFilter(datr, norder)/(norder^2);
    
    datr = preprocess.removeChannelArtifacts(datr, ops.artifactThresh, ops.artifactNchans, 50);

    datr    = datr * Wrot; % whiten the data and scale by 200 for int16 range

    if debug
        d = datr;
        subplot(212)
        imagesc(d')
        xlim([0 1e3])
        colormap jet
        drawnow
%         plot(d + (1:size(d,2))*1300, 'k')
    end
    
    
    datcpu  = gather(int16(datr)); % convert to int16, and gather on the CPU side
    fwrite(fidW, datcpu, 'int16'); % write this batch to binary file
end

rez.Wrot    = gather(Wrot); % gather the whitening matrix as a CPU variable

fclose(fidW); % close the files
fclose(fid);

fprintf('Time %3.0fs. Finished preprocessing %d batches. \n', toc, Nbatch);

rez.temp.Nbatch = Nbatch;
