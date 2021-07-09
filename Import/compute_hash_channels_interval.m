function [hash,smallfp] = compute_hash_channels_interval(DataFolder,ShankMap,Interval)
%*******
%** function hash = compute_hash_channels_interval(DataFolder,ShankMap,Interval)
%***
%*** input:
%***     DataFolder:  folder containing INTAN continuous files
%***     ShankMap:  if not [], then map for channels 
%***     Interval:  if [], not used, else start and end time in Ephys secs
%***
%*** output:
%***   hash struct
%***   always uses common mode references if more than a single channel
%***    contains sp:  per channel, hash spike times
%***   smallfp struct
%***     contains LFP but only for epochs specified by interval (smaller)
%**************
  
  tic
  %**********
  hash.DataFolder = DataFolder;
  hash.sp = [];
  hash.ChanNum = 0;
  hash.Map = [];
  hash.BaseRate = 20;
  %******* quick LFP over one interval
  smallfp.DataFolder = DataFolder;
  smallfp.LINE_NOISE_FREQ = [60 120 180];
  smallfp.LINE_NOISE_WID = 1.0;
  smallfp.Time = [];
  smallfp.LFP = [];
  
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
  smallfp.LFP = cell(1,ChanNum);
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
      
       %********** select subset of relevant time for filtering
       if (k == 1)
          if ~isempty(Interval)
             tims = find( (sptime >= Interval(1)) & (sptime < Interval(2)) );                 
          else
             tims = 1:length(sptime);
          end
          if ~isempty(tims)
            %****** select sampling and filters
            Fs = spinfo.header.sampleRate;
            Ksamp = floor( Fs / 1000);  % spacing for ~1Khz sampling
            samples = 1:Ksamp:length(spdata(tims));
            [lowb,lowa] = butter(1,(300/Fs*2),'low');
            %******* filtering for noise (on down-sampled data)
            NoFs = (Fs/Ksamp);
            bandb = cell(1,length(smallfp.LINE_NOISE_FREQ));
            banda = cell(1,length(smallfp.LINE_NOISE_FREQ));
            for zk = 1:length(smallfp.LINE_NOISE_FREQ)
              [bandb{zk},banda{zk}] = butter(3, [((smallfp.LINE_NOISE_FREQ(zk)-1)/NoFs*2), ...
                                           ((smallfp.LINE_NOISE_FREQ(zk)+1)/NoFs*2)]);
            end
          end
       end
       %**********
       if isempty(tims)
           continue;
       else
           disp(sprintf('Filtering subset %d of %d total',length(tims),length(spdata)));
       end
       %*******
       disp(sprintf('Low-passing for LFP from continuous file %s',tag));
       hzdata = filter(lowb,lowa,spdata(tims));
       hzdata = flipud(hzdata);
       hzdata = filter(lowb,lowa,hzdata);
       hzdata = flipud(hzdata);
       %******
       disp('Downsampling data');
       hzdata = hzdata(samples);
       %******
       disp('Noise filtering');
       ndata = hzdata;
       for zk = 1:length(smallfp.LINE_NOISE_FREQ)
          hzdata = filter(bandb{zk},banda{zk},hzdata);
          hzdata = flipud(hzdata);
          hzdata = filter(bandb{zk},banda{zk},hzdata);
          hzdata = flipud(hzdata);
          ndata = ndata - hzdata;  % remove line noise
       end
       %******* storage of time and LFP
       if (k == 1)
         smallfp.Time = sptime(tims(samples));
       end
       smallfp.LFP{kmap} = ndata;
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
  savname = [DataFolder,filesep,tagname,'_hashint'];
  disp(sprintf('Saving hash struct to %s',savname));
  save(savname,'hash','-v7.3');
  disp('Hash struct saved');
  savname = [DataFolder,filesep,tagname,'_smallfp'];
  disp(sprintf('Saving hash struct to %s',savname));
  save(savname,'smallfp','-v7.3');
  disp('Struct is saved, processing done!');
  disp(' ');
  toc
  
return;