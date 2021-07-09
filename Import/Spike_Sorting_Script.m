%%********* SPIKE SORTING SCRIPT ********

%% ****
CN = ChanNums;
if ischar(ShankLayout)
    MAP = load(ShankLayout);
else
    MAP = 1:ChanNums;
end

Clust = num2cell(MAP);

%% Spike Thresholding on selected sub-channels
%%********* data loaded, now you would do spike thresholding on it
%********* to make a list of spike times
sp = cell(1,CN);
for k = 1:CN
    %***** first average the channels from cluster
    disp('Spike sortings on channels:');

    disp('Averaging channels ...');
    for j = 1:size(Clust{k},2)
       tag = sprintf('*CH%d.continuous',Clust{k}(j));
       SpikeFiles = dir([DataFolder,filesep,tag]);
       if isempty(SpikeFiles)
           fprintf('Error finding continuous file %s\n',tag);
           continue;
       end 
       [spdata,sptime,spinfo] = read_ephys.load_open_ephys_data_faster([DataFolder,filesep,SpikeFiles(1).name]);
       if (j == 1)
           avsp = spdata;
       else
           avsp = avsp + spdata;
       end
    end
    %********************
    sp{k} = spikesort.single_channel_threshold(avsp,spinfo, 1);
    sp{k}.st = sptime( sp{k}.ss );
    %*********************
end
%****************************************************************

