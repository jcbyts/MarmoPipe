function igood = get_dead_channels(ops, chanMap)
% find dead channels by looking for correlations between neighboring channels
% that are way off

Nbatch = ops.Nbatch;
twind = ops.twind;
NchanTOT = ops.NchanTOT;
NT = ops.NT;
Nchan = numel(chanMap);
load(ops.chanMap);

fid = fopen(ops.fbinary, 'r');
% irange = [NT/8:(NT-NT/8)];

ibatch = 1;

v = zeros(Nchan, 1);

% 2D differentiating filter
k = [0 -1 0; -1 4 -1; 0 -1 0];
CC = 0; % initialize 

% from a subset of batches, count threshold crossings
while ibatch<=Nbatch
    
    offset = twind + 2*NchanTOT*NT* (ibatch-1);
    fseek(fid, offset, 'bof');
    buff = fread(fid, [NchanTOT NT], '*int16');

    if isempty(buff)
        break;
    end
    
    data = single(buff)';
    data = data - median(data,2);
    CC = CC + cov(data,1);
    
    v = v + var(single(buff(chanMap,:)), [], 2)/Nbatch;
    ibatch = ibatch + ops.nSkipCov;
end

CC = CC / ceil((Nbatch-1)/ops.nSkipCov);

sqd = conv2(CC, k, 'same').^2;


A = (xcoords(chanMap)-xcoords(chanMap)')+100;
sqmask = conv2(A,k,'same');
sqd = sqd .* (sqmask.^2==0);
sqd = sqd .* (1 - eye(size(sqd)));
imagesc(sqd)

chbadness = mean(sqd,2)' + mean(sqd);

figure(1); clf
plot(chbadness); hold on
plot(xlim, median(chbadness)*[1 1]*10, '--r')
drawnow

igood = true(Nchan,1);
thresh = median(v)*[2 1/2];
igood(v < thresh(2) | v > thresh(1)) = false;
igood(chbadness > 10*median(chbadness)) = false;

fclose(fid);
figure, plot(v)
hold on
plot(find(~igood), v(~igood), 'o')
plot(xlim, thresh(1)*[1 1], '--')
plot(xlim, thresh(2)*[1 1], '--')
title('Dead Channels')
xlabel('Channel')
ylabel('Variance')

fprintf('found %d bad channels \n', sum(~igood))
drawnow