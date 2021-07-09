%*********** ANALYSIS SCRIPT *******************************
%***********************************************************
%*** NOTE:  The goal here is just to analyze a single file
%***        This will not be very well organized
%***        Once we get an initial analysis done, we will
%***        build a GUI to do this systematically and
%***        provide a larger structure
%***       
%***        For now, all info is imported to an Exp struct
%***        which is saved back to the data folder and then
%***        can be loaded for any analysis
%***
%***        Revised so it concatenates all files from a session
%***
%***        Things still to do:
%***        At present, we are still not:
%***          1) storing information about electrode type
%***          2) storing information about LFP (notch filtered)
%***          3) using Kilosort for multi-channel data sorting
%***          4) storing session info to Meta-Table
%************************************************************
%***
%***  1/2/19:   Modified Import_Script_DDPI.m to include code to 
%***            read .edf files if they are present, Jude
%***

%% Spike File importing

marmoPath = addMarmoPipe();

user = 'gravedigger';
SERVER_DATA_DIR = addKiloPaths(user);
%% Select folder and convert raw data to a single binary file
if ~exist('DataFolder', 'var')
    assert(exist('SERVER_DATA_DIR', 'var')==1, 'Kilosort2ImportScript: SERVER_DATA_DIR does not exist. you must point to the raw data before running this')
    DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
end

% if contains(DataFolder,'Logan')
%     
%     sName = ('NandyFrontFlip');
%     %shank = hardware.electrodeFactory('NandyFrontFlip');%<-- CHECK IF THIS IS RIGHT
%     shank = hardware.electrodeFactory(sName);%<-- CHECK IF THIS IS RIGHT
%     disp('Logan selected')
%     disp(sName)
% else
shank = hardware.electrodeFactory('Nandy64');
%     shank = hardware.electrodeFactory('Nandyflip');
%     shank = hardware.electrodeFactory('ShankA'); % 32 channel probe
% end


io.oe2dat(DataFolder, shank, 'overwrite', false);
ops = io.loadOps(DataFolder);
SavePath = ops.root; % location of the processed data
currDir = pwd;

%% 
%************
TAGSTART = 63;   
TAGEND = 62; 
%*************
user = 'jakegravedigger';

switch user
    case 'judehome'
        SERVER_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Raw';
        PROCESSED_DATA_DIR = 'C:\Users\jmitchell\Dropbox\FovealTrans\DataAnalysis\DDPI_Processed';
    case 'shannawork'
        SERVER_DATA_DIR = 'C:\Users\Shanna\Dropbox\Marmo Lab Website\PSA\DDPI_Raw';
        PROCESSED_DATA_DIR = 'C:\Users\Shanna\Dropbox\Marmo Lab Website\PSA\DDPI_Processed';
    case 'jakegravedigger'
%         SERVER_DATA_DIR = 'C:\PSA\Raw3'; %'C:\Raw';
%         PROCESSED_DATA_DIR = 'C:\PSA\Processed3';%'Z:\Data\Processed'; %'C:\Processed';
        SERVER_DATA_DIR = 'C:\Raw';
        PROCESSED_DATA_DIR = 'C:\PSA\Processed3';%'Z:\Data\Processed'; %'C:\Processed';
end

% DataFolder = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
[~, FileTag] = fileparts(DataFolder);
 DataFolder = fullfile(SERVER_DATA_DIR,FileTag);





%% plot the raw just to check that the channel map is correct
inds = (1:30e3) + 55000e3+55000e3; % sample index

inVolts = true;
data = io.loadRaw(ops, inds, inVolts);
data = bsxfun(@minus, data, mean(data)); % common average reference
figure(1); clf
channelOffsets =(1:ops.Nchan)*200;
plot(bsxfun(@plus, data', channelOffsets))
% title(SavePath)

%% run Kilosort without GUI
mitchell_kilosort
 
%% Run manual merge stage in phy
fprintf('Once you finish sorting in the GUI and you have saved the output\n')
fprintf('Open Anaconda Prompt and enter\n')
fprintf('activate phy2\n')
fprintf('cd %s\n', SavePath)
fprintf('phy template-gui params.py\n')

%% save spikes as struct
%cd(currDir);
ops = io.loadOps(DataFolder);
info = io.loadEphysInfo(DataFolder);
sp = io.getSpikesFromKilo(ops, info);
% sp = loadKSdir(SavePath);
save(fullfile(SavePath, 'spikes-kilo.mat'), '-v7.3', '-struct', 'sp')


%% Get spikes from Kilosort
prl = dir(fullfile(DataFolder, '_processed*'));
prl(arrayfun(@(x) ~x.isdir, prl)) = []; % remove any non-directories
npr = numel(prl);
fprintf('Found %d processed folders\n', npr);
sp = {};
ctr = 1;
for i = 1:npr
    spkFiles = dir(fullfile(DataFolder, prl(i).name, '*spikes*.mat'));
    if ~isempty(spkFiles)
        sp = load(fullfile(DataFolder, prl(i).name, spkFiles(1).name));
    end
end

%**********************
osp = sp;  % save old
cids = unique(sp.clu);
N = numel(cids);
sp = cell(1, N);
for k = 1:N
    sp{k}.st = osp.st(osp.clu==cids(k));
end

%% **** based on what you set, it runs spike sorting (eventually KiloSort)
% search data directory for sorted spike files, and use those if available
if (1)  % skip this assuming you did Kilosort
 spfiles = dir([SERVER_DATA_DIR,filesep,FileTag,filesep,'*.*_spk.mat']);
 if (~isempty(spfiles))
    fprintf('Using .spk files present in data directory for spike times\n');
    disp('Counting flagged units ....');
    CN = 0;
    holdsp = cell(1,1);
    for zk = 1:size(spfiles,1)
       fname = [SERVER_DATA_DIR,filesep,FileTag,filesep,spfiles(zk).name];
       load(fname);
       if ~isempty(sp)
         for k = 1:size(sp,2)
             CN = CN + 1;
             holdsp{1,CN} = sp{1,k};
             holdsp{1,CN}.name = [fname,'_',sprintf('U%d',k)];
         end
       end
    end
    clear sp;
    sp = holdsp;
    fprintf('Total of %d units identified in sorted files\n',CN);
  else
    % ******** Call the Spike Sorting Script with parameters set *******
    SingleChannel = false;  % if false the specify shank layout
    ChNumber = NaN;   % if single channel, set here the number
    ShankLayout = 'ShankA_32map.txt';  % shank layout name
    ChanNums = 32;  % if a shank recording, how many channels, 1 to N
    %*******
    Spike_Sorting_Script;   % otherwise use a simple thresholding scrip
 end
end
%** returns sp ... a struct with spike times and other info

%% Events File importing
%******** now grab the events file with strobes
EventFiles = dir([DataFolder,filesep,'*.events']);
if isempty(EventFiles)
    disp('Error finding events file');
    return;
end
[evdata,evtime,evinfo] = read_ephys.load_open_ephys_data([DataFolder,filesep,EventFiles(1).name]);
%**** convert events into strobes with times
[tstrobes,strobes] = read_ephys.convert_data_to_strobes(evdata,evtime,evinfo);
strobes = read_ephys.fix_missing_strobes(strobes);
disp('Stobes are loaded');


%% Loading up the MarmoView Data Files
%****************************************************************
ExpFiles = dir([DataFolder,filesep,'*z.mat']);
if isempty(ExpFiles)
    disp('Error finding *z.mat file');
    return;
end
%****** get order by date of files to import them
FileDates = cellfun(@datenum,{ExpFiles(:).date});
DoFileSort = [ (1:size(FileDates,2))' FileDates'];
FileSort = sortrows(DoFileSort,2); % sort them by date
%***** read and append files in order by date *********
BigN = 0;
for zk = FileSort(:,1)'
  fname = ExpFiles(zk).name;  
  load([DataFolder,filesep,fname]);
  if ~BigN
    Exp.D = D;
    Exp.S = S;
    BigN = size(D,1); % number of trials
  else
    for k = 1:size(D,1)
       Exp.D{BigN+k} = D{k};  % appending trials 
    end  
    BigN = BigN + size(D,1);
  end
  clear D;
  clear S;
  fprintf('Experiment file %s loaded\n',fname);
end
%***** store spikes info in Exp struct, let's keep all info there
%***** once finished we will clear everything but the Exp struct
if exist('osp', 'var')
    Exp.osp = osp;
end
Exp.sp = sp;
% for k = 1:size(sp,2)  % now done inside spike sorting
%   Exp.sp{k}.st = sptime( Exp.sp{k}.ss );
% end
clear sp;
disp('Experiment files loaded');
%**************************

%% Synching up strobes from Ephys to MarmoView
%******************* returns start and end times in ephys record, or NaN if missing
%***** this code would look simple, but sometimes one bit is flipped in
%***** error and then you have to play catch up to find the missing code
%***** since there is some redundancy (start and end codes for each trial)
%***** this gives us a way to recover cases with just one errant bit
disp('Synching up ephys strobes');
for k = 1:size(Exp.D,1)
   start = synchtime.find_strobe_time(TAGSTART,Exp.D{k}.STARTCLOCK,strobes,tstrobes);
   finish = synchtime.find_strobe_time(TAGEND,Exp.D{k}.ENDCLOCK,strobes,tstrobes);
   if (isnan(start) || isnan(finish))
       fprintf('Synching trial %d\n',k);
       if (isnan(start) && isnan(finish))  % if both are missing drop the trial
             fprintf('Dropping entire trial %d from protocol %s\n',k,Exp.D{k}.PR.name);
             Exp.D{k}.START_EPHYS = NaN;
             Exp.D{k}.END_EPHYS = NaN;          
       else
           %******* here we could try to repair a missing code, or find a
           %******* a partial code nearby
           Exp.D{k}.START_EPHYS = start;
           Exp.D{k}.END_EPHYS = finish;
           tdiff = Exp.D{k}.eyeData(end,6) - Exp.D{k}.eyeData(1,1);
           if isnan(start) && ~isnan(finish)
               disp('**Approximating start code from end');
               Exp.D{k}.START_EPHYS = Exp.D{k}.END_EPHYS - tdiff;
               %****** now see if you can do even better
               %****** see if the real code is there but a bit flipped
               zz = find(tstrobes == finish);  % find end code
               istart = zz(1) - 7;
               if (istart >= 1) && (strobes(istart) == TAGSTART)  % candidate start before end
                 beftag = strobes((istart+1):(istart+6))';
                 mato = sum( Exp.D{k}.STARTCLOCK & beftag );             
                 if (mato >= 5)  % all but one of taglet matched
                   E.D{k}.START_EPHYS = tstrobes(istart);
                   disp('****Located matching start strobe, one bit was flipped');
                 end
               end
               %*******************************
           end
           if isnan(finish) && ~isnan(start)
               disp('##Approximating end code from start');
               Exp.D{k}.END_EPHYS = Exp.D{k}.START_EPHYS + tdiff;
               %****** now see if you can do even better
               %****** see if the real code is there but a bit flipped
               zz = find(tstrobes == start);  % find end code
               istart = zz(1) + 7;
               if (istart < size(strobes,1)) && (strobes(istart) == TAGEND)  % candidate start before end
                 endtag = strobes((istart+1):(istart+6))';
                 mato = sum( Exp.D{k}.ENDCLOCK & endtag );
                 if (mato >= 5)  % all but one of taglet matched
                   E.D{k}.END_EPHYS = tstrobes(istart);
                   disp(endtag)
                   disp('####Located matching end strobe, one bit was flipped');
                 end
               end
               %*******************************
           end
           %****************
       end
   else
      Exp.D{k}.START_EPHYS = start;  % otherwise store the NaN -- maybe trial not used
      Exp.D{k}.END_EPHYS = finish;
   end
end
disp('Finished Synching up ephys strobes');

%% Loading up the VPX file as a long data stream
% Can be done later, might just use matlab eye data for now
%*****************************************************************
VPX = 0;
DDPI = 0; 
EDF = 0;
VpxFiles = dir([DataFolder,filesep,'*.vpx']);
if isempty(VpxFiles)
    % try to find ddpi files if no VPX files
    VpxFiles = dir([DataFolder,filesep,'*.ddpi']); 
    if isempty(VpxFiles)
       %******** if not DDPI, then try to find .edf from EyeLink  
       VpxFiles = dir([DataFolder,filesep,'*edf.mat']);
       if (isempty(VpxFiles))
           disp('Error finding raw eye data file');
           return;
       else
           EDF = 1;
       end
    else
        DDPI = 1;
    end
else
    VPX = 1;
end
%****** get order by date of files to import them
FileDates = cellfun(@datenum,{VpxFiles(:).date});
DoFileSort = [ (1:size(FileDates,2))' FileDates'];
FileSort = sortrows(DoFileSort,2); % sort them by date
%****************************
if (DDPI == 1)  % compute median eye params
   %****** to screen bad ddpi samples, feed in median values
   %****** for raw eye val to eye dva, screen any values
   %****** outside +/- 20 degs and faster than 2dva steps
   pixPerDeg = Exp.S.pixPerDeg;
   cc = [];
   dxx = [];
   dyy = [];
   for ii = 1:numel(Exp.D)
          cc = [cc ; Exp.D{ii}.c];
          dxx = [dxx ; Exp.D{ii}.dx];
          dyy = [dyy ; Exp.D{ii}.dy];
   end
   c = median(cc);
   dx = median(dxx);
   dy = median(dyy);
   %***************
end
%***** read and append files in order by date *********
BigN = 0;
for zk = FileSort(:,1)'
  fname = VpxFiles(zk).name;  
  vpx_filename = [DataFolder,filesep,VpxFiles(zk).name];
  if (DDPI == 1)
      vpx = read_ddpi.load_ddpi_file(vpx_filename);%,pixPerDeg,c,dx,dy);  % makes DDPI look like VPX
  else
      if (EDF == 1)
          vpx = read_edf.load_edf_file(vpx_filename);
      else
          vpx = read_vpx.load_vpx_file(vpx_filename);
      end
  end
  if ~isempty(vpx)
      if ~BigN
        Exp.vpx = vpx;
      else
        %******* append time to prevent overlap across files  
        vpx.raw(:,1) = vpx.raw(:,1) + BigN;  % add time offset before concat
        vpx.smo(:,1) = vpx.smo(:,1) + BigN;
        if isfield(vpx,'dat')
           vpx.dat(:,1) = vpx.dat(:,1) + BigN;
        end
        vpx.tstrobes = vpx.tstrobes + BigN;
        %******* concatenate large file stream **********
        Exp.vpx.raw = [Exp.vpx.raw ; vpx.raw];
        Exp.vpx.smo = [Exp.vpx.smo ; vpx.smo];
        if isfield(vpx,'dat')
           Exp.vpx.dat = [Exp.vpx.dat ; vpx.dat];
        end
        Exp.vpx.tstrobes = [Exp.vpx.tstrobes ; vpx.tstrobes];
        Exp.vpx.strobes = [Exp.vpx.strobes ; vpx.strobes];
        %******** compute new last time, plus one minute
      end
      BigN = Exp.vpx.smo(end,1) + 60.0; % last time point plus one minute
      clear vpx;
      disp('******************************************');
      disp(sprintf('Experiment file %s loaded',fname));
  else
      disp(sprintf('WARNING: failed to read %s',fname));
  end
end


%% Synching up the strobes from VPX to MarmoView
% Same thing, might use matlab eye data for now
%*****************************************************************
disp('Synching up vpx strobes');
for k = 1:size(Exp.D,1)
   start = synchtime.find_strobe_time(TAGSTART,Exp.D{k}.STARTCLOCK,Exp.vpx.strobes,Exp.vpx.tstrobes);
   finish = synchtime.find_strobe_time(TAGEND,Exp.D{k}.ENDCLOCK,Exp.vpx.strobes,Exp.vpx.tstrobes);
   if (isnan(start) || isnan(finish))
       disp(sprintf('Synching VPX trial %d',k));
       if isnan(finish) && isnan(start)
           disp(sprintf('Dropping entire VPX trial %d from protocol %s',k,Exp.D{k}.PR.name));
           Exp.D{k}.START_VPX = NaN;
           Exp.D{k}.END_VPX = NaN;
       else
           %******* here we could try to repair a missing code, or find a
           %******* a partial code nearby
           Exp.D{k}.START_VPX = start;
           Exp.D{k}.END_VPX = finish;
           tdiff = Exp.D{k}.eyeData(end,6) - Exp.D{k}.eyeData(1,1);
           if isnan(start) && ~isnan(finish)
               Exp.D{k}.START_VPX = Exp.D{k}.END_VPX - tdiff;
               disp('Approximating VPX start code');
           end
           if isnan(finish) && ~isnan(start)
               Exp.D{k}.END_VPX = Exp.D{k}.START_VPX + tdiff;
               disp('Approximating VPX end code');
           end
           %****************
       end
   else
      Exp.D{k}.START_VPX = start;
      Exp.D{k}.END_VPX = finish;
   end
end
disp('Finished synching up vpx strobes');

%% Saccade processing:
%** Perform basic processing of eye movements and saccades
%** and link up smoothed eye position with list of saccade times
Saccade_Script;
disp('Saccade trials are done processing');


%% Save phy clusters (good, mua, noise)
Exp.PHYClu = readtable([ops.root '\cluster_group.tsv'], 'FileType', 'text');

%%  compute the CSD onsets if available
CSD_Spike_Onset_Script;  % determine onsets of CSD task
minCSD = min(CSD.Onsets_Ephys)-0.5;  % first CSD flash in secs
maxCSD = max(CSD.Offsets_Ephys)+0.5; % last CSD flash in sec
Exp.CSD = CSD;
 
%% Save Exp Struct to data folder for future use (PROC folder)
%****** should consider later how we want to store things (per protocol,
%****** but off hand it seems feasible given the generic structure of
%****** of the D struct maybe we could concatenate all sessions in one
ExpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'a.mat'];

Exp.ProcDataFolder = PROCESSED_DATA_DIR;
Exp.DataFolder = DataFolder;
Exp.FileTag = FileTag;
save(ExpFile,'-v7.3','Exp');
fprintf('Exp struct saved to %s\n',ExpFile);
fprintf('Mat File %s\n',FileTag);


%% HASH sort all the channels and store a struct
[hash,smallfp] = compute_hash_channels_interval(DataFolder,[],[minCSD,maxCSD]);
hashFile = [PROCESSED_DATA_DIR,filesep,FileTag,'_hashint.mat'];
save(hashFile,'-v7.3','hash');
fprintf('Hash struct saved to %s\n',hashFile);
smallfpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'_smallfp.mat'];
save(smallfpFile,'-v7.3','smallfp');
fprintf('Smallfp struct saved to %s\n',smallfpFile);
fprintf('Mat File %s\n',[FileTag,'_smallfp']);
disp('Saving Hash and Smallfp done');


%******* FINAL PROCESSING OF LFP AND HASH IS SLOW ... go get food etc...
%% LFP
%**** run with true to create the LFP data in folders
%**** note, only a partial struct is returned
%[lfp,~,lfpinfo] = io.getLFP_cmode(ops,false,false);
[lfp, ~, lfpinfo] = io.getLFP(ops, true, false);
lfpinfo.ProcDataFolder = PROCESSED_DATA_DIR;
lfpinfo.DataFolder = DataFolder;
lfpinfo.FileTag = FileTag;
lfpFile = [PROCESSED_DATA_DIR,filesep,FileTag,'_lfp.mat'];
save(lfpFile,'-v7.3','lfp','lfpinfo');
fprintf('LFP struct saved to %s\n',lfpFile);
fprintf('Mat File %s\n',[FileTag,'_lfp']);
disp('Saving LFP done');

%% at the end clear it up 
%*********** clear up environment when finished  
clear all;
  
