function Wrot = get_whitening_matrix(rez)
% based on a subset of the data, compute a channel whitening matrix
% this requires temporal filtering first (gpufilter)

ops = rez.ops;
Nbatch = ops.Nbatch;
twind = ops.twind;
NchanTOT = ops.NchanTOT;
NT = ops.NT;
NTbuff = ops.NTbuff;
chanMap = ops.chanMap;
Nchan = rez.ops.Nchan;
xc = rez.xc;
yc = rez.yc;

% load data into patches, filter, compute covariance
if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
    [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
end

% -------------------------
% USE JRCLUST FILTERING
filtOpts.freqLimBP = [ops.fshigh ops.fslow];
filtOpts.freqLimNotch = [];
filtOpts.freqLimStop = [];
filtOpts.nSamplesPad = 100;
filtOpts.useGPUFilt = true;
filtOpts.sampleRate = ops.fs;
filtOpts.gainBoost = 1;
filtOpts.useElliptic = true;
filtOpts.filtOrder = 3;
norder = 2;
nchannel = 32;
nskip = 2;
nav = 3;
Icar = toeplitz([zeros(1,nskip) ones(1,nav) zeros(1,nchannel-nav-nskip)]);
if ops.Nchan > 32
    Icar = blkdiag(Icar,Icar);
end
Icar = Icar(chanMap,chanMap);
Icar = Icar ./ sum(Icar);


fprintf('Getting channel whitening matrix... \n');
fid = fopen(ops.fbinary, 'r');
CC = gpuArray.zeros( Nchan,  Nchan, 'single'); % we'll estimate the covariance from data batches, then add to this variable

NS = 0;
ibatch = 1;
debug = getOr(ops, 'debug', 0);



while ibatch<=Nbatch
    offset = max(0, twind + 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));
    fseek(fid, offset, 'bof');
    buff = fread(fid, [NchanTOT NTbuff], '*int16');

    if isempty(buff)
        break;
    end
    nsampcurr = size(buff,2);
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
    end

%     datr    = gpufilter(buff, ops, rez.ops.chanMap); % apply filters and median subtraction
    
    

    datr = filtfiltChain(single(buff(chanMap,:))', filtOpts);
    datr = doFFTClean(datr, 5, false);
    datr = bsxfun(@minus, datr, cast(single(datr)*Icar, 'like', datr) ); % common average reference

    datr = jrclust.filters.sgFilter(datr, norder)/(norder^2);
        
    siteRMS = jrclust.utils.estimateRMS(datr, 10e5);
    qqFactor = 4;
    siteThresh = siteRMS*qqFactor;

    [i,~] = find(datr < -siteThresh);
    nt = size(datr, 1);
    nospks = setdiff(1:nt, unique(i + (-30:30)));
    
    CC        = CC + (datr(nospks,:)' * datr(nospks,:))/NT; % sample covariance
%     CC        = CC + (datr' * datr)/NT; % sample covariance

    ibatch = ibatch + ops.nSkipCov; % skip this many batches
end
CC = CC / ceil((Nbatch-1)/ops.nSkipCov); % normalize by number of batches

fclose(fid);

if ops.whiteningRange<Inf
    if ops.whiteningRange < 0 % whitening by shank (defined along xc)
        eps = 1;
        xcs = unique(xc);
        nx = numel(xcs);
        ccs = cell(nx);
        for ixc = 1:nx
            ccs{ixc} = whiteningFromCovariance(CC(xc==xcs(ixc), xc==xcs(ixc)), eps);
        end
        Wrot = blkdiag(ccs{:});
    else
        % if there are too many channels, a finite whiteningRange is more robust to noise in the estimation of the covariance
        ops.whiteningRange = min(ops.whiteningRange, Nchan);
        Wrot = whiteningLocal(gather(CC), yc, xc, ops.whiteningRange); % this function performs the same matrix inversions as below, just on subsets of channels around each channel
    end
else
    Wrot = whiteningFromCovariance(CC, 1/NT);
end

Wrot    = ops.scaleproc * Wrot; % scale this from unit variance to int 16 range. The default value of 200 should be fine in most (all?) situations.

fprintf('Channel-whitening matrix computed. \n');
