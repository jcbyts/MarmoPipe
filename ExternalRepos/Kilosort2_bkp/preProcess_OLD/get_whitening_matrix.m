function Wrot = get_whitening_matrix(rez)

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

fprintf('Getting channel whitening matrix... \n');
fid = fopen(ops.fbinary, 'r');
if ops.GPU
    CC = gpuArray.zeros( Nchan,  Nchan, 'single');
else
    CC = zeros( Nchan,  Nchan, 'single');
end


% irange = [NT/8:(NT-NT/8)];

ibatch = 1;
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
    if ops.GPU
        dataRAW = gpuArray(buff);
    else
        dataRAW = buff;
    end
    dataRAW = dataRAW';
    dataRAW = single(dataRAW);
    dataRAW = dataRAW(:, chanMap);
    dataRAW = dataRAW(:,ops.igood);
    
    % subtract the mean from each channel
    dataRAW = dataRAW - mean(dataRAW, 1);
    
    % CAR, common average referencing by median
    if getOr(ops, 'CAR', 1)
        dataRAW = dataRAW - median(dataRAW, 2);
    end
    
%     datr = fft(dataRAW, [], 1);
%     datr = datr./(1 + abs(datr));
%     datr(irange, :) = 0;
%     datr = real(ifft(datr, [], 1));
%     
    datr = filter(b1, a1, dataRAW);
    datr = flipud(datr);
    datr = filter(b1, a1, datr);
    datr = flipud(datr);

     [datr, bad] = preprocess.removeChannelArtifacts(datr, ops.artifactThresh, ops.artifactNchans, 50);
%     % CAR, common average referencing by median
%     if getOr(ops, 'CAR', 1)
%         datr = datr - median(datr, 2);
%     end
    
    nospikes = zscore(sum(datr.^2,2))<15;
    nospikes(bad) = false;
    CC(ops.igood,ops.igood)        = CC(ops.igood,ops.igood) + (datr(nospikes,:)' * datr(nospikes,:))/NT;    
    
    ibatch = ibatch + ops.nSkipCov;
end
CC = CC / ceil((Nbatch-1)/ops.nSkipCov);
CC = CC + eye(size(CC,1));
fclose(fid);
fprintf('Channel-whitening filters computed. \n');

if ops.whiteningRange<Inf
    ops.whiteningRange = min(ops.whiteningRange, Nchan);
    Wrot_ = whiteningLocal(gather_try(CC(ops.igood,ops.igood)), yc(ops.igood), xc(ops.igood), ops.whiteningRange);
    Wrot = zeros(size(CC), 'like', CC);
    Wrot(ops.igood,ops.igood) = Wrot_;
else
    [E, D] 	= svd(CC(ops.igood,ops.igood));
    D       = diag(D);
    eps 	= mean(D); %1e-6;
%     eps 	= 1e-6;
    f  = mean((D+eps) ./ (D+1e-6));
    fprintf('%2.2f ', f)
%     Wrot_ 	= E * diag(f./(D + eps).^.5) * E';
    Wrot_ 	= E * diag(f./(D + eps)) * E';
    Wrot = zeros(size(CC), 'like', CC);
    Wrot(ops.igood,ops.igood) = Wrot_;
    
end
Wrot    = ops.scaleproc * Wrot;
