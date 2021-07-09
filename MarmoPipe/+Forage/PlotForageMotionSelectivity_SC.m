function PlotForageMotionSelectivity(Exp,SPClust)
% function PlotForageMotionSelectivity(Exp,SPClust)
%   input: takes Exp struct as the input
%  Exp fields:
%    D: {100x1 cell}  - data per trial
%    S: [1x1 struct]  - rig settings struct
%   sp: [1x1 struct]  - spikes from spike sorting
%   cp: (will be a 60hz filtered LFP at 1000 hz)
%   vpx: a struct for eye position
%  SPClust:  which unit to use in analysis
%
%   output: plot of spike raster sorted by motion direction
 
  Smoothing = 10.0; % Gauss smoothing for PSTH plot
  %******** Parameters for windows on spike counts per trial
  BefStim = 100;  % these are in ms, but below seconds are the spk and eye times
  AftStim = 400;  % so will convert to ms at some point
  PostStim_Start = 50;  % interval to compute counts from, start
  PostStim_End =  350;  % and end in ms
  
  
  %***** first identify which trials had motion stims, and directions
  Tlist = [];
  MoList = [];
  for zk = 1:size(Exp.D,1)
     if (length(Exp.D{zk}.PR.name) >= 6)
        if strcmp(Exp.D{zk}.PR.name(1:6),'Forage')
          if (Exp.D{zk}.PR.noisetype == 4) % Motion
              Tlist = [Tlist ; zk];
              nh = Exp.D{zk}.PR.NoiseHistory;
              zz = find( nh(:,2) > 0 ); % stim onset frames
              MoList = [MoList ; nh(zz,2)];  % get histogram of dirs sampled
              disp(sprintf('Registering trial %d',zk));
              
              %******* COMMENTS ON FORMAT
              % NoiseHistory: column 1 is the time (matlab)
              %               column 2 is zero, or positive for stim onset
              %                          and if an onset then an integer
              %                          to represent the motion number
              %                        or negative, for stim offset with
              %                          a corresponding negative motion
              %                          number
              %   we will aim to time lock on motion onsets for now
          end
        end
     end
  end
  MoNum = max(MoList);
  vh = 1:MoNum;
  MoFreq = hist(MoList,vh);
  
  %**** start some plotting
  FH = figure;
  subplot('position',[0.6 0.6 0.3 0.3]);
  MoDirs = (vh - 1) * 360 / MoNum;
  bar(MoDirs,MoFreq,'b'); hold on;
  bmax = max(MoFreq);
  axis([-15 360 0 (bmax * 1.2)]);
  xlabel('Motion Number');
  ylabel('Trial Repeats');
  
  %**************
  sp = Exp.sp{SPClust};
  %*********************************************
  
  %******** loop through all trials and build a raster plot time-locked
  %******** on the stimulus onset
  Rast = cell(1,MoNum);
  RastCnt = zeros(1,MoNum);
  RastList = cell(1,MoNum);
  for k = 1:MoNum   % similar to PlotActivityAroundSaccade in Flag
                    % but here doing the raster per motion condition
     Rast{k} = zeros(MoFreq(k),(BefStim+AftStim+1));
  end
  %**** trial by trial stats and counts
  PostDir = [];
  PostSpk = [];
  %***********
  for zk = 1:size(Tlist,1)   % loop through relevant trials
      %*******
      tr = Tlist(zk);
      matstart = Exp.D{tr}.eyeData(1,1);
      %***
      nh = Exp.D{tr}.PR.NoiseHistory;
      zz = find( nh(:,2) > 0 ); % stim onset frames
      mdirs = nh(zz,2);  % motion directions for onsets
      mtims = nh(zz,1);  % frame times in matlab total time
      mtims = mtims - matstart; % time relative to trial start
      %******* find relevant times in EPHYS record around each event
      for ek = 1:size(mdirs,1)   % loop over events
        %***********
        modir = mdirs(ek);
        RastCnt(modir) = RastCnt(modir) + 1;
        mocnt = RastCnt(modir);
        %*********** now find spikes in epoch
        staspk = Exp.D{tr}.START_EPHYS + mtims(ek) - (BefStim/1000); % covert to secs
        finspk = Exp.D{tr}.START_EPHYS + mtims(ek) + (AftStim/1000);
        z = find( (sp.st > staspk) & (sp.st < finspk) );
        %***********
        if ~isempty(z)  % load up spikes to rasters
          difftime = ceil( 1000 * (sp.st(z) - staspk ) );  % integers in ms
          Rast{modir}(mocnt,difftime) = 1;
          RastList{modir} = [RastList{modir} ; [(difftime-BefStim) ones(size(difftime))*mocnt]];
        end
        %**** get spike counts per event and store them (fastest plot)
        %*********** now find spikes in post-stim epoch
        post_staspk = Exp.D{tr}.START_EPHYS + mtims(ek) - (PostStim_Start/1000); % covert to secs
        post_finspk = Exp.D{tr}.START_EPHYS + mtims(ek) + (PostStim_End/1000);
        z = find( (sp.st > post_staspk) & (sp.st < post_finspk) );
        %***********
        if ~isempty(z)
            rate = length(z) * (1000/(PostStim_End - PostStim_Start));
        else
            rate = 0;
        end
        %******
        PostDir = [PostDir ; modir];
        PostSpk = [PostSpk ; rate];
        %******
      end
      %**********
      disp(sprintf('Processing trial %d of %d',tr,size(Tlist,1)));
      %*******************
  end
  
  %**** plot the spike counts per motion direction (simple report to start)
  figure(FH);
  subplot('position',[0.6 0.1 0.3 0.3]);
  PostDirs = (PostDir - 1) * 360 / MoNum;
  SpkMax = max(PostSpk);
  plot(PostDirs,PostSpk,'ko');
  axis([-15 360 0 (SpkMax * 1.2)]);
  xlabel('Motion Direction');
  ylabel('Spike Counts');
  
  
return;

function smo = gauss_smooth(psth,Gsig)

    % Make the number of samples depending on the gaussian window size
    gaussian_filter_size = 2*Gsig-1; % if Gsig = 10, 19 samples total
                                     % 9 left & 9 right from the mean

    % Make smoothing kernel using gaussian filter
    for i = 1:gaussian_filter_size
        gauss  = exp(-(((i-Gsig).^2)./(2*Gsig^2)));
        gauss_filter(i,:) = gauss;
    end
    % Normalize the gaussian filter
    gauss_smooth = gauss_filter/sum(gauss_filter);

    psth_size    = length(psth);
    filter_size  = length(gauss_smooth);
    filter_cent = floor((filter_size+1)/2);

    for i=1:psth_size   % size_smooth

        % Always 0 for the initial value (only sum from product of two vectors)
        smo(i) = 0;
        nomo(i) = 0;

        % Apply filter to data
        for j = 1:filter_size
             diff = (j-filter_cent);   % this coordinate, goes - 3*sig up to 3*sig
             samp = (i+diff);
             if ( (samp >= 1) && (samp <= psth_size) )
                 smo(i) = smo(i) + (psth(samp) * gauss_smooth(j));
                 nomo(i) = nomo(i) + gauss_smooth(j);
             end       
        end
        %********
        if (nomo(i) > 0)
            smo(i) = smo(i) / nomo(i);
        end
        %***********
    end
return;

function [u_smooth,sem_smooth] = compute_psth(Raster,Smoothing)
% function [u_smooth,sem_smooth] = compute_psth(Raster,Smoothing)
%   computes the mean and Jacknife error bars
%  input:
%     Raster is Trials x Time points in ms
%     Smoothing is Gaussian sigma in ms
   
   %******* compute the Jacknife
   smoothsub = [];
   N = size(Raster,1);
   JNum = 20;   % the number of Jacknife's in estimate
   for i = 1:JNum 
    aex = 1+floor((i-1)*N/JNum);
    bex = ceil(i*N/JNum);
    excludeset = [1:aex,bex:N];
    psth = mean(Raster(excludeset,:));
    psth = psth*1000;
    smooth_data = gauss_smooth(psth,Smoothing);  
    smoothsub = [smoothsub; smooth_data];
   end
   u_smooth = mean(smoothsub);
   sem_smooth = (JNum-1) * var(smoothsub);
   sem_smooth = sqrt(sem_smooth);   % variance is multiplied by N-1 Jacknife

return;
