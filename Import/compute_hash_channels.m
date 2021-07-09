function hash = compute_hash_channels(DataFolder,ShankMap)
%*******
%** function hash = compute_hash_channels(DataFolder,ShankMap)
%***
%*** input:
%***     DataFolder:  folder containing INTAN continuous files
%***     ShankMap:  if not [], then map for channels 
%***
%*** output:
%***   hash struct
%***   always uses common mode references if more than a single channel
%***    contains sp:  per channel, hash spike times
%***    
%**************
  
  tic
  %**********
  hash.DataFolder = DataFolder;
  hash.sp = [];
  hash.ChanNum = 0;
  hash.Map = [];
  hash.BaseRate = 20;
  
  %********* find all continuous files
  spfiles = dir( [DataFolder,filesep,'*CH*continuous*'] );
  ChanNum = length(spfiles);
  if (ChanNum == 0)
      return;
  else
      hash.ChanNum = ChanNum;
  end
  
  %******* check if ShankMap matches the number of files
  if ~isempty(ShankMap)
      MapNum = length(ShankMap);
      if (MapNum < ChanNum)
          disp('Fewer shank maps than number of continuous files');
          hash.Map = [];
      else
          hash.Map = ShankMap;
      end
  else
      hash.Map = [];
  end
  
  %****** if more than one channel, compute average over channels 
  %****** (so to use common mode referencing)
  AvChan = [];
  for k = 1:ChanNum
       tag = [DataFolder,filesep,spfiles(k).name];
       disp(sprintf('Computing ch %d common mode average from continuous file %s',k,tag));
       [spdata,sptime,spinfo] = read_ephys.load_open_ephys_data_faster(tag);
       if (k == 1)
           AvChan = spdata;
       else
           AvChan = AvChan + spdata;
       end
  end
  if ~isempty(AvChan)
      AvChan = AvChan / ChanNum;   % compute average
  end
  disp('Average common mode computed');

  hash.sp = cell(1,ChanNum);
  %****** read in each channel, common mode ref, and filter for MUA, LFP,
  %****** and for hash spikes
  for k = 1:ChanNum
       fname = spfiles(k).name;
       tag = [DataFolder,filesep,fname];
       %*** search for channel in tag name
       kchan = 0;
       for ii = 2:length(fname)
           if strcmp('CH',fname((ii-1):ii))
               jj = ii+1;
               while (~isempty(str2num(fname(jj))))
                   kchan = str2num(fname(jj)) + (kchan*10);
                   jj = jj + 1;
               end
           end
       end
       if (kchan == 0)
           disp(sprintf('Unable to find channel in name %s',fname));
           continue;
       else
           disp(sprintf('\n Processing channel %d from name %s',kchan,fname));
       end
       %*******
       kmap = kchan;  %you can implement remapping of channels later
       [spdata,sptime,spinfo] = read_ephys.load_open_ephys_data_faster([DataFolder,filesep,fname]);
       %*** common mode referenced if you have it
       if ~isempty(AvChan)
           spdata = spdata - AvChan;
       end
       
       %********* now filter the raw data **************
       disp(sprintf('Hash spike sorting continuous file %s for %d',tag,kmap));
       hash.sp{kmap} = spikesort.auto_channel_threshold(spdata,spinfo, hash.BaseRate);
       hash.sp{kmap}.st = sptime( hash.sp{kmap}.ss );      
      %***************************************************
  end
  toc
  %****** save result for storage, do not return it (too large?)
  tic
  ksep = 0;
  LK = length(DataFolder);
  for k = 1:LK
      if (DataFolder(LK+1-k) == filesep)
         ksep = (LK+1-k);
         break;
      end
  end
  if (ksep)
      tagname = DataFolder((ksep+1):LK);
  else
      tagname = 'noname';
  end
  savname = [DataFolder,filesep,tagname,'_hash'];
  disp(sprintf('Saving hash struct to %s',savname));
  save(savname,'hash','-v7.3');
  disp('Struct is saved, processing done!');
  disp(' ');
  toc
  
return;