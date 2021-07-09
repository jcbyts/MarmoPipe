%*** just practicing with Gabe
%*** run through trials in Exp events struct
%*** find those trials with CSD and then
%*** find the onset times, and time lock spiking raster
%*** to those times to validate we are using time correctly
%***
%*** Gabe will work on integrating onset times into CSD analysis then

sanitycheck = 0;  % if you want to see debugging

%*** Find all the CSD trials in the Exp struct
Trials = length(Exp.D);
CSD_Trials = [];
for k = 1:Trials
   ProtoName = Exp.D{k}.PR.name;
%    if (strcmp(ProtoName,'ForageProceduralNoise'))  % for Jake
   if (strcmp(ProtoName,'Forage'))   % for Shanna
       NoiseType = Exp.D{k}.PR.noisetype;
       disp(sprintf('Trial(%d); %s  NoiseType(%d)',k,ProtoName,NoiseType));
       if (NoiseType == 6)
           CSD_Trials = [CSD_Trials ; k];
       end
   end
end

%%
%*** now loop through CSD trials and get all the onset times
%*** of stimuli, and put them in spike reference time
CSD_Onsets = [];
CSD_Offsets = [];
CSD_Onsets_Ephys = [];
CSD_Offsets_Ephys = [];
NTrials = length(CSD_Trials);
for k = 1:NTrials
    kk = CSD_Trials(k);
    NoHist = Exp.D{kk}.PR.NoiseHistory;
%     CSD_Onsets2 = [CSD_Onsets2; Exp.D{kk}.START_EPHYS];
    %*** Noise History is time in column 1, and contrast (0 or 1) in col 2
    %** find all the Onsets, as transition from 0 to 1 in column 2
    for i = 2:size(NoHist,1)
       if (NoHist(i-1,2) == 0) && (NoHist(i,2) >= 1)   % 
           tt = NoHist(i,1);
           CSD_Onsets = [CSD_Onsets ; tt];  % store Mat time
           %******* convert to Ephys time per trial start-clocks
%            tt = tt - Exp.D{kk}.STARTCLOCKTIME;  % Jake - 0 from start of trial
          tt = tt - Exp.D{kk}.eyeData(6,1);    % Shanna - start of trial in mat time
           tt = tt + Exp.D{kk}.START_EPHYS;     % time from start of trial in ephys
           CSD_Onsets_Ephys = [CSD_Onsets_Ephys ; tt];
           % search for corresponding offset, if trial ends, insert a NaN
           for j = (i+1):size(NoHist,1)
             if (NoHist(j-1,2) >= 1) && (NoHist(j,2) == 0)   % 
               tt = NoHist(j,1);
               CSD_Offsets = [CSD_Offsets ; tt];  % store Mat time
               %******* convert to Ephys time per trial start-clocks
%                tt = tt - Exp.D{kk}.STARTCLOCKTIME;  % 0 from start of trial
             tt = tt - Exp.D{kk}.eyeData(6,1);    % Shanna - start of trial in mat time 
               tt = tt + Exp.D{kk}.START_EPHYS;     % time from start of trial in ephys
               CSD_Offsets_Ephys = [CSD_Offsets_Ephys ; tt];
               break;
             end
           end
           if (j == size(NoHist,1))  % Offset occured after trial end
               CSD_Offsets = [CSD_Offsets ; NaN];
               CSD_Offsets_Ephys = [CSD_Offsets_Ephys ; NaN];
           end
       end 
    end
    %*** sanity check result looks right per trial
    if (sanitycheck == 1)
        figure(10); hold off;
        plot(NoHist(:,1),NoHist(:,2),'k-'); hold on;
        plot(CSD_Onsets,(0.5*ones(size(CSD_Onsets))),'rx');
        xlabel('Time (secs)');
        input('check');
    end
end

%%
%***** what you need for aligning LFP is CSD_Onsets_Ephys
%****** but lets make a spike raster to confirm we got the timing
%****** correct, because there should be a visual response to it
if (1) %(sanitycheck == 1)
 for Unit = 1:size(Exp.sp,2) 
   NOnsets = length(CSD_Onsets_Ephys);
   SpkRaster = [];
   OffRaster = [];
   SpChan = Unit;
   for k = 1:NOnsets
     tt = CSD_Onsets_Ephys(k);
     z = find( (Exp.sp{SpChan}.st >= tt) & (Exp.sp{SpChan}.st < (tt+0.5) ) );
     if ~isempty(z)
      sptt = Exp.sp{SpChan}.st(z) - tt;  % time lock spikes relative to onset
      SpkRaster = [SpkRaster ; [(k*ones(size(sptt))) sptt]];
     end
     ott = CSD_Offsets_Ephys(k) - tt;
     OffRaster = [OffRaster ; [k ott]];
   end
   if ~isempty( SpkRaster )
     figure(10); hold off;
     plot(1000*SpkRaster(:,2),SpkRaster(:,1),'k.'); hold on;
     plot(1000*OffRaster(:,2),OffRaster(:,1),'r.'); hold on;
     xlabel('Time (ms)');
     ylabel('Trials');
     title(sprintf('CSD Raster: Unit(%d)',SpChan));
     input('check');
   else
     figure(10); hold off;
     
     disp(sprintf('Not plotting for unit %d, no spikes',Unit));
    end
 end
end

