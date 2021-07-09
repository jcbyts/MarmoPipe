

%% LFP Processing script
%% 0
addMarmoPipe()

%% 1

pick_session = false;

%************
TAGSTART = 63;   
TAGEND = 62;  
%*************
user = 'bluethunder';

switch user
    case 'judehome'
        SERVER_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Raw';
        PROCESSED_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Processed';
    case 'bluethunder'
%         SERVER_DATA_DIR = 'C:\PSA_Gravedigger\Raw';
        SERVER_DATA_DIR = 'C:\Raw';
%         PROCESSED_DATA_DIR = 'C:\PSA_Gravedigger\Processed';
        PROCESSED_DATA_DIR = 'Z:\Data\Processed_Laminar'; %C:\Processed';
    case 'jakegravedigger'
        SERVER_DATA_DIR = 'C:\Raw';
        PROCESSED_DATA_DIR = 'Z:\PSA_EDF\EDF_Processed_Decoding\X-Y Recordings';
end

if pick_session
    DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
else 
%     Sessions = {'Milo_2021-02-15_12-41-31_MT32_1b', 'Milo_2021-02-22_11-26-14_MT32_4b', 'Milo_2021-03-17_11-57-28_MT64_5b', 'Milo_2021-03-18_11-05-37_MT64_6b', 'Milo_2021-04-07_12-40-34_MT64_7b', 'Milo_2021-04-12_11-54-37_MT64_9c', 'Milo_2021-04-16_12-15-05_MT64_11b', 'Milo_2021-04-21_11-51-59_MT64_14b', 'Milo_2021-05-10_11-49-28_MT64_20b', 'Milo_2021-05-14_11-02-57_MT64_22b', 'Milo_2021-05-26_13-01-25_MT64_26b', 'Milo_2021-06-08_11-19-32_MT64_29b'};
    Sessions = {'Milo_2021-03-17_11-57-28_MT64_5b', 'Milo_2021-03-18_11-05-37_MT64_6b', 'Milo_2021-04-07_12-40-34_MT64_7b', 'Milo_2021-04-16_12-15-05_MT64_11b', 'Milo_2021-04-21_11-51-59_MT64_14b', 'Milo_2021-05-10_11-49-28_MT64_20b', 'Milo_2021-05-26_13-01-25_MT64_26b', 'Milo_2021-06-08_11-19-32_MT64_29b'};
    sess_num = 1;
    sess = Sessions{sess_num};
    disp(['Session: ' sess])
    DataFolder = sess;
end
[~, FileTag] = fileparts(DataFolder);
sess = FileTag;
% DataFolder = fullfile(SERVER_DATA_DIR,FileTag);


if ~exist(char(strcat(PROCESSED_DATA_DIR,filesep,FileTag, '.mat')), 'file')
    tag_split = split(FileTag, '_');
    subj_name = tag_split(1);
    date = split(tag_split(2), '-');
    yr = split(date(1), '0');
    yr = yr(2);
    day = date(3);
    month = date(2);
    FileTag = strcat(subj_name, '_', day, month, yr);
    FileTag = FileTag{1};
    disp(['FileTag: ', FileTag])
end

if ~exist(char(strcat(PROCESSED_DATA_DIR,filesep,FileTag, '.mat')), 'file')
    error('experiment file not found')
end

disp('done')




%% 2 THIS IS A HACK TO CHANGE OPS PATH - MAY NOT NEED IN PIPELINE
ops = io.loadOps(DataFolder);
% [fPath, fName, fExt] = fileparts(ops.root);
% newStr = split(ops.root,'/');
newStr=regexp(ops.root,'\','split');
new_root = fullfile(SERVER_DATA_DIR,newStr{end-1},newStr{end});
ops.root = new_root;
disp(new_root)

newStr=regexp(ops.fbinary,'\','split');
new_fbinary = fullfile(new_root,newStr{end});
ops.fbinary = new_fbinary;
disp(new_fbinary)

%% 3 Load EXP file
disp('loading exp file')
ExpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'.mat'];
load(ExpFile);
disp('done')

%% 4 LFP

% xcoords = Exp.osp.xcoords;
% ycoords = Exp.osp.ycoords;
%**** run with true to create the LFP data in folders
%te, only a partial struct is returned
lfpFile = fullfile(PROCESSED_DATA_DIR,'lfp',[FileTag '_lfp.mat']);

if ~exist(lfpFile, 'file') || pick_session
    [data, timestamps, info] = io.getLFP(ops, true, false);
else
    warning('File already exists.. loading existing file. ')
    warning('Delete old file to rerun LFP processing. ')
    disp('loading lfp struct')
    lfp = load(lfpFile);
    data = lfp.data;
    timestamps = lfp.timestamps;
    info = lfp.info;
    disp('done')
end
unique_x = unique(Exp.osp.xcoords);
num_shanks = size(unique_x, 1);
shank_len = size(data, 2)/num_shanks;
deadChan = []; % Add dead channels here if there are any (1-32 first shank, 33-64 second shank, etc.)
if num_shanks > 0
    
    xcoords = ones(shank_len,num_shanks);
    for i = 1:num_shanks
        xcoords(:,i) = xcoords(:,i)*unique_x(i);
    end
end
diff_ycoords = abs(mode(diff(Exp.osp.ycoords)));
ycoords = repmat(flip(linspace(diff_ycoords, diff_ycoords*shank_len, shank_len)), num_shanks, 1)';

%lfpFile = fullfile(PROCESSED_DATA_DIR,'lfp',[FileTag '_lfp.mat']);
% lfpFile = [PROCESSED_DATA_DIR,filesep,'lfp',filesep,FileTag,'_lfp.mat'];
disp('Saving LFP struct');
save(lfpFile, '-v7.3', 'timestamps', 'info', 'data', 'xcoords', 'ycoords', 'deadChan')
fprintf('LFP struct saved to %s\n',lfpFile);
fprintf('Mat File %s\n',[FileTag,'_lfp']);
disp('Saving LFP done');

%% 5 Import LFP (once saved)
disp('loading lfp struct')
lfpFile = fullfile(PROCESSED_DATA_DIR,'lfp',[FileTag '_lfp.mat']);
lfp = load(lfpFile);
disp('done')

%% 6 PLOT CSD - noisetype 3 (full field flash)
stats = csd.getCSD(lfp, Exp, 'plotIt', false, 'noisetype', 3);
figure(1); clf
csd.plotCSD(stats, 'overlayLFP', false)
title([FileTag '_noistype3'],'Interpreter','none')
xlabel('Time (ms)')
ylabel('Depth')


%% 6 PLOT CSD - noisetype 6 (moving flashes)
stats = csd.getCSD(lfp, Exp, 'plotIt', false, 'noisetype', 6);
figure(2); clf
csd.plotCSD(stats, 'overlayLFP', false)
title([FileTag '_noistype6'],'Interpreter','none')
xlabel('Time (ms)')
ylabel('Depth')

%% at the end clear it up
%*********** clear up environment when finished
clear all;
close all;

