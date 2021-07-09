
%%

% FileTag = 'Ellie_2019-09-17_13-53-24_MT64D9';
% 
% SERVER_DATA_DIR = 'C:\Raw';
% PROCESSED_DATA_DIR = 'C:\Processed';
% 
% % DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
% DataFolder = [SERVER_DATA_DIR,filesep,FileTag];

Nchan = 64;
nSamples = 10e3;
data = zeros(Nchan, nSamples);

for i = 1:Nchan
    fprintf('Loading channel %d\n', i)    
    chFiles = dir(fullfile(DataFolder, sprintf('*CH%d.continuous', i)));
    if ~isempty(chFiles)
        filename = fullfile(DataFolder, chFiles.name);
        
        data_seg = load_open_ephys_data_faster(filename, 'unscaledInt16');
        data(i,:) = double(data_seg(5500e3 + (1:nSamples)));
    end
end
    
%%

data = dataorig;
%%


% shank = hardware.electrodeFactory('Nandy');
shank{1} = hardware.electrode.NeuroNexus64_2_32;
shank{1}.headstages{1} = hardware.headstage.intan_RHD2164_2;
% shank{1}.headstages{2} = hardware.headstage.intan_RHD2132;
        
% data = data - median(data);
headstage = hardware.combineHeadstages(shank{1}.headstages);

channelMap = headstage.mapChannels(shank{1});


sprintf('%d,', channelMap)
% channelMap(channelMap>32) = channelMap(channelMap>32)
% channelMap = headstage.channelMap(1:64);
% channelMap = 1:Nchan;

figure(1); clf
channelOffsets =(1:Nchan)*1000;
% plot(bsxfun(@plus, data(channelMap,:)', channelOffsets))

imagesc(data(channelMap,:))
% plot(bsxfun(@plus, double(squeeze(DATA(:,:,1))), 10*channelOffsets(1:61)))
% [map.chanMap(:) map.connected(:)]
d = data(channelMap,:)';
CC = corr(d - median(d,2));
figure(2);clf, subplot(121), imagesc(cov(data')), subplot(122), imagesc(CC)

% CC = CC );
load(ops.chanMap)
%%
k = [0 -1 0; -1 4 -1; 0 -1 0];
sqd = conv2(CC, k, 'same').^4;
figure(1); clf
A = (xcoords(channelMap)-xcoords(channelMap)')+100;
sqmask = conv2(A,k,'same');
sqd = sqd .* (sqmask.^2==0);
sqd = sqd .* (1 - eye(size(sqd)));
imagesc(sqd, .5*[-1 1])

chbadness = mean(sqd,2)' + mean(sqd);

plot(chbadness); hold on
plot(xlim, median(chbadness)*[1 1]*10, '--r')

inds = find(chbadness > 10*median(chbadness));
clf
% plot(CC(inds,:)')
plot(CC(:,inds))

d = data(channelMap,:);
% clf
% plot(d(setdiff(1:64, inds),:)', 'Color', .4*[1 1 1]); hold on
% plot(d(inds,:)')

% plot(mean(diff(sqd),2) + mean(diff(sqd,[],2),1)')% + mean(sqd))
%%
[u, s, v] = svd(CC);

n = 2;
xdiff = (CC - u(:,1:n)*s(1:n,1:n)*v(:,1:n)').^2;

imagesc(xdiff)
%%
[~, ind] = sort(mean(xdiff.* (1- eye(size(CC)))), 'descend');
plot(mean(xdiff.* (1- eye(size(CC)))))
ind(1:3)
%%
d1 = data(channelMap(1:32),:);

figure, imagesc(corr(d1'))
%%
[a,b]=sortBatches2(1-corr(d1'));
figure(1); clf
imagesc(a)
%%
%% save channel map to file
% offset channels > 32
channelMap(channelMap > 32) = channelMap(channelMap > 32) + 3;
%%


for i = 1:numel(channelMap)
    fprintf('%d) %d\n', i, channelMap(i))
end

