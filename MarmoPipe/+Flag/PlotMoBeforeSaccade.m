function PlotMoActivityAroundSaccade(Exp,SPClust)
% function PlotActivityAroundSaccade(Exp)
%   input: takes Exp struct as the input
%  Exp fields:
%    D: {100x1 cell}  - data per trial
%    S: [1x1 struct]  - rig settings struct
%   sp: [1x1 struct]  - spikes from spike sorting
%   cp: (will be a 60hz filtered LFP at 1000 hz)
%   vpx: a struct for eye position
%
%   output: plot of spike raster aligned to saccades to targets in Camoflag
 
  %*** it will go through all Camoflag trials and determine 
  %*** the saccades into the aperture
  EyeRadius = 2.5;  % should be defined in protocol, but what we count
                    % as the saccade end point landing close enough to 
                    % the target
  AmpThresh = 0.6;  % fraction of saccade distance of total target ecc
                    % before we include that saccade as totally correct
  Smoothing = 10.0; % Gauss smoothing for PSTH plot
  TooFast = 0.10;    % throw out saccades with RT too fast
  %***********************************
  
  %******** Parameters for windows on spike counts per trial
  BefSac = 150;  % these are in ms, but below seconds are the spk and eye times
  AftSac = 300;  % so will convert to ms at some point
  %****** selects windows for analysis (non-overlapping now)
  PreWin = [-20,40];
  StimWin = [80,160];
  %**********
  XJit = 15; %7.5;  % orientation jit
  YJit = 0.1; % percentage rate
  %********************************
  %**************
  sp = Exp.sp{SPClust};
  
  %****** experiment specific saccade flagging (if needed)
  SacNum = 0;
  CorNum = 0;
  FixNum = 0;
  Tlist = [];  % list of trials with saccades correct
  TargList = [];  % list of which target he made a saccade towards
  TPosList = [];  % list of target choice, and x,y of sac endpoint
  for i = 1:size(Exp.D,1)
        if  strcmp(Exp.D{i}.PR.name(1:6),'FlagMo')
         %if (1) % (Exp.D{i}.PR.deltaOri == 0)  % isnan(Exp.D{i}.PR.changori)  
         if ~isfield(Exp.D{i}.PR,'target_null') || (Exp.D{i}.PR.target_null ~= 1)  % no blank target in RF location
          if (Exp.D{i}.PR.error == 0) || (Exp.D{i}.PR.error == 4)  % only correct saccade trials (for now) 
               eyeSmo = Exp.D{i}.eyeSmo;
               slist = Exp.D{i}.slist;
               FixRadius = Exp.D{i}.P.fixWinRadius;  % within this bound for be fixating 
               
               %***********
               Tx = Exp.D{i}.PR.targ_x;
               Ty = Exp.D{i}.PR.targ_y;
               TargDist = norm([Tx(1),Ty(1)]);  % eccentricity of target
               BestTarg = NaN;
               BestEtx = NaN;
               BestEty = NaN;
               
               if ~isempty(slist)
                  foundfix = 0;
                  foundsac = 0;
                  foundcor = 0;
                  for k = 1:size(slist,1)
                      stasac = slist(k,4); % integer start time of saccade
                      endsac = slist(k,5); % integer end time of saccade
                      stx = eyeSmo(stasac,2);
                      sty = eyeSmo(stasac,3);
                      etx = eyeSmo(endsac,2); 
                      ety = eyeSmo(endsac,3);
                      fdist = norm([etx,ety]);
                      if (fdist < FixRadius)   % was animal in fixation
                          foundfix = 1;
                      end
                      if (foundfix == 1)
                         edist = [];
                         for ik = 1:size(Tx,1)
                            edist(ik) = norm([Tx(ik),Ty(ik)]-[etx,ety]);
                         end
                         if (min(edist) < EyeRadius)  % first saccade into aperture
                            slist(k,7) = 1;  % mark the saccade
                            foundcor = 1;
                            %****** compute fraction distance of saccade
                            sacamp = norm([(etx-stx),(ety-sty)]);
                            if (sacamp > (TargDist * AmpThresh) )
                               slist(k,7) = 2; % flag the saccade was large
                               foundsac = 1;
                               zz = find( edist == min(edist));
                               BestTarg = zz(1);   % which was target
                               BestEtx = etx;
                               BestEty = ety;
                               break;
                            end
                         end
                      end
                  end
                  %**************
                  offstim = Exp.D{i}.PR.stimOffset - Exp.D{i}.eyeData(1,1);
                  z = find( Exp.D{i}.eyeData(:,5) == 4);
                  if ~isempty(z)
                    onstim = Exp.D{i}.eyeData(z(1),1) - Exp.D{i}.eyeData(1,1);
                  else
                    onstim = NaN;
                  end
                  myoffstim = onstim + 0.1;
                  if ( abs(offstim-myoffstim) > 0.01)
                      offstim = myoffstim;
                  end
                  %*******
                  z = find( Exp.D{i}.eyeData(:,5) == 5);  % first moment Matlab detects leaving fix
                  if ~isempty(z)
                     msactime = Exp.D{i}.eyeData(z(1),1) - Exp.D{i}.eyeData(1,1);
                  else
                     disp(sprintf('No saccade onset found trial %d',i)); 
                     msactime = NaN;
                  end
                  %***************
                  offsac = find( slist(:,7) == 2);  % find flagged saccade
                  if ~isempty(offsac)
                      ontime = slist(offsac(1),1);  % onset of that saccade
                      if ( (offstim-ontime) < 0.03)
                         if ~isnan(msactime) && (msactime > ontime)     
                             FixNum = FixNum + foundfix;
                             CorNum = CorNum + foundcor;         
                             SacNum = SacNum + foundsac;
                             if (foundsac == 1) 
                                 if ~isnan(msactime) && ~isnan(onstim)
                                   if ( (msactime-onstim) > TooFast)
                                      Tlist = [Tlist ; i];
                                      TargList = [TargList ; BestTarg];
                                      %*********************
                                      TPosList = [TPosList ; [BestTarg,BestEtx,BestEty]];
                                      %************
                                   end
                                 end
                             end
                         end
                      end
                  end
                  %********
               end
               %********* store new saccade list
               Exp.D{i}.slist = slist;
               %*******************
           end
          end
        end
  end
  disp(sprintf('Fixations %d, Correct %d, Saccades Flagged %d',FixNum,CorNum,SacNum)); 
  %*********************************************
  disp(sprintf('Presaccadic, %d att and %d unatt',size(find(TargList == 1),1),...
                                                  size(find(TargList ~= 1),1)));
  
  %******** loop through all trials and build a raster plot time-locked
  %******** on the saccade offset (should see some visual response??)
  StimRast = zeros(SacNum,(BefSac+AftSac+1));  % all zeros and ones
  StimRastList = [];
  StimSacDur = [];
  StimMSacDur = [];
  %*******************
  SacRast = zeros(SacNum,(BefSac+AftSac+1));  % all zeros and ones
  SacRastList = [];
  SacSacDur = [];
  SacMSacDur = [];
  %*******************
  %**** trial by trial stats and counts
  zPreOri = [];
  PreSpk = cell(1,1);
  StimSpk = cell(1,1); 
  %**********
  tr = Tlist(1);
  ChangList = [];
  %***********
  for k = 1:size(Tlist,1) % SacNum
     tr = Tlist(k);  % trial number
     tg = TargList(k);
     %*******
     slist = Exp.D{tr}.slist;
     offsac = find( slist(:,7) == 2);  % find flagged saccade
     ontime = slist(offsac(1),1);  % onset of that saccade
     offtime = slist(offsac(1),2);  % saccade offset in seconds from trial start
     %****** find the stimulus onset, which is first step into state 2
     z = find( Exp.D{tr}.eyeData(:,5) == 4);
     if ~isempty(z)
       onstim = Exp.D{tr}.eyeData(z(1),1) - Exp.D{tr}.eyeData(1,1);
     else
       onstim = NaN;
     end
     %******** stim offset time
     offstim = Exp.D{tr}.PR.stimOffset - Exp.D{tr}.eyeData(1,1);
     %******* not working well, off-stim
     StimSacDur = [StimSacDur ; [((1000*(ontime-onstim))) k]]; 
     SacSacDur = [SacSacDur ; [((1000*(ontime-ontime))) k]];
     %****** compute stim onset *****
     z = find( Exp.D{tr}.eyeData(:,5) == 5);  % Matlab detects leaving fix
     if ~isempty(z)
        msacdur = Exp.D{tr}.eyeData(z(1),1) - Exp.D{tr}.eyeData(1,1);
     else
        msacdur = NaN; 
     end
     %****
     if ~isnan(msacdur)
        StimMSacDur = [StimMSacDur ; [((1000*(msacdur-onstim))) k]];
        SacMSacDur = [SacMSacDur ; [((1000*(msacdur-ontime))) k]];
     else
        StimMSacDur = [StimMSacDur ; [NaN NaN]];
        SacMSacDur = [SacMSacDur ; [NaN NaN]];
     end
    
     %****************
     for zk = 1:2
         if (zk == 1)
             TimeLockOff = 2;  % lock from stim onset
         else
             TimeLockOff = 0;  % lock from sac onset
         end

         if (TimeLockOff == 2)
            staspk = Exp.D{tr}.START_EPHYS + onstim - (BefSac/1000); % covert to secs
            finspk = Exp.D{tr}.START_EPHYS + onstim + (AftSac/1000);  
         else
            staspk = Exp.D{tr}.START_EPHYS + ontime - (BefSac/1000); % covert to secs
            finspk = Exp.D{tr}.START_EPHYS + ontime + (AftSac/1000);  
         end
         %**** get spike times from the interval of analysis
         z = find( (sp.st > staspk) & (sp.st < finspk) );
         if ~isempty(z)
             difftime = ceil( 1000 * (sp.st(z) - staspk ) );  % integers in ms
             if (zk == 1)
                 StimRast(k,difftime) = 1;
                 StimRastList = [StimRastList ; [(difftime-BefSac) ones(size(difftime))*k]];
             else
                 SacRast(k,difftime) = 1;
                 SacRastList = [SacRastList ; [(difftime-BefSac) ones(size(difftime))*k]];                    
             end
         end

         if (TimeLockOff == 2)
              PreSta = StimWin(1);
              PreEn = StimWin(2);
              onstaspk = Exp.D{tr}.START_EPHYS + onstim + (PreSta/1000); % covert to secs
              onfinspk = Exp.D{tr}.START_EPHYS + onstim + (PreEn/1000);
              z = find( (sp.st >= onstaspk) & (sp.st < onfinspk) );
              if ~isempty(z)
                 StimSpk{1} = [StimSpk{1} ; (length(z)/(onfinspk-onstaspk)) ];
              else
                 StimSpk{1} = [StimSpk{1} ; 0];
              end
         else
              PreSta = PreWin(1);
              PreEn = PreWin(2);
              onstaspk = Exp.D{tr}.START_EPHYS + ontime + (PreSta/1000); % covert to secs
              onfinspk = Exp.D{tr}.START_EPHYS + ontime + (PreEn/1000);
              z = find( (sp.st >= onstaspk) & (sp.st < onfinspk) );
              if ~isempty(z)
                 PreSpk{1} = [PreSpk{1} ; (length(z)/(onfinspk-onstaspk)) ];
              else
                 PreSpk{1} = [PreSpk{1} ; 0];
              end
         end
     end
     %**********
     disp(sprintf('Processing trial %d of %d',k,SacNum));
     %******* this gets more complicated now, as depends on which target
     if ~isnan(tg) 
          ori = Exp.D{tr}.PR.targ_motion(tg);
          dori = ori + Exp.D{tr}.PR.deltaOri;
          if (dori > 360)
             dori = dori - 360;
          end
          if (dori < 0)
             dori = dori + 360;
          end 
          preori = Exp.D{tr}.PR.targori;  % always in RF
          apreori = ori;
          postori = dori;  % change trials 
     else
          preori = NaN;
          apreori = NaN;
          postori = NaN;
     end
     %*************
     if (tg == 1)  % target in RF
         ChangList = [ChangList ; 1]; % att trial
     else
         ChangList = [ChangList ; 2];
     end
     %******** get spike counts pre and post -saccade
     zPreOri =  [zPreOri  ; preori ];  % ori in RF bef sac
     %*******************
  end
  %******** now figure out what is the unique set of orientations
  iz = ~isnan(zPreOri);
  Orivals = unique(zPreOri(iz));
  OriNum = length(Orivals);
  
  %*********
  PreOri = [];
  for k = 1:size(zPreOri,1)
      pr_ori = find( zPreOri(k) == Orivals ); 
      if ~isempty(pr_ori)
         PreOri =  [PreOri  ; pr_ori(1)];  % assumes 12 orientations
      else
         PreOri = [PreOri ; NaN];
      end
  end
  %***********
  TrA = find(ChangList == 1);
  TrU = find(ChangList == 2);
  disp(sprintf('Into RF %d and Away %d',size(TrA,1),size(TrU,1)));
  %****************
  
  %****** plot the stim raster ****************
  H = figure;
  set(H,'Position',[100 100 1200 800]);
  subplot('position',[0.05 0.4 0.25 0.5]);
  %********
  zb = find( ChangList == 1);
  zbb = ismember(StimRastList(:,2),zb);
  hh = plot(StimRastList(zbb,1),StimRastList(zbb,2),'r.'); hold on;
  set(hh,'Markersize',0.2);
  %******
  zc = find( ChangList == 2);  % blank trials
  zcc = ismember(StimRastList(:,2),zc);
  hh = plot(StimRastList(zcc,1),StimRastList(zcc,2),'b.'); hold on;
  set(hh,'Markersize',0.2);
  %******
  hh = plot(StimSacDur(:,1),StimSacDur(:,2),'k.'); 
  hh2 = plot(StimMSacDur(:,1),StimMSacDur(:,2),'g.');
  %*******
  axis([-BefSac AftSac 0 SacNum]);
  plot([0,0],[0,SacNum],'b-');
  xlabel('Time (ms)');
  ylabel('Trials');
  title('(Timelock Stim Onset)'); 
  %***************************************
  [Auu,Asu] = compute_psth(StimRast(TrA,:),Smoothing);
  [Uuu,Usu] = compute_psth(StimRast(TrU,:),Smoothing); 
  %***** plot the psth *******************
  hold on;
  subplot('position',[0.05 0.05 0.25 0.25])
  %****** subset of trials where no ori change occured
  h2 = plot(-BefSac:AftSac,Auu,'r-');hold on;
  h2b = plot(-BefSac:AftSac,Auu+(2*Asu),'r:');hold on;
  h2b = plot(-BefSac:AftSac,Auu-(2*Asu),'r:');hold on;
  %****** subset of trials where ori change DID occur
  h2 = plot(-BefSac:AftSac,Uuu,'b-');hold on;
  h2b = plot(-BefSac:AftSac,Uuu+(2*Usu),'b:');hold on;
  h2b = plot(-BefSac:AftSac,Uuu-(2*Usu),'b:');hold on;
  %************
  axis tight;
  V = axis;
  axis([-BefSac AftSac 0 V(4)]);
  plot([0,0],[0,V(4)],'k-');
  plot([StimWin(1),StimWin(1)],[0,V(4)],'m-');
  plot([StimWin(2),StimWin(2)],[0,V(4)],'m-');
  xlabel('Time (ms)');
  ylabel('Rate (hz)');
  title('PSTH');
  
  
  subplot('position',[0.4 0.4 0.25 0.5]);
  %********
  zb = find( ChangList == 1);
  zbb = ismember(SacRastList(:,2),zb);
  hh = plot(SacRastList(zbb,1),SacRastList(zbb,2),'r.'); hold on;
  set(hh,'Markersize',0.2);
  %******
  zc = find( ChangList == 2);  % blank trials
  zcc = ismember(SacRastList(:,2),zc);
  hh = plot(SacRastList(zcc,1),SacRastList(zcc,2),'b.'); hold on;
  set(hh,'Markersize',0.2);
  %******
  hh = plot(SacSacDur(:,1),SacSacDur(:,2),'k.');      % leave in vpx
  hh2 = plot(SacMSacDur(:,1),SacMSacDur(:,2),'g.');   % leave fix in matlab
  %***********
  axis([-BefSac AftSac 0 SacNum]);
  plot([0,0],[0,SacNum],'b-');
  xlabel('Time (ms)');
  ylabel('Trials');
  title('(Timelock Onset)');   
  %***************************************
  [Auu,Asu] = compute_psth(SacRast(TrA,:),Smoothing);
  [Uuu,Usu] = compute_psth(SacRast(TrU,:),Smoothing);
  %***** plot the psth *******************
  hold on;
  subplot('position',[0.4 0.05 0.25 0.25])
  %****** subset of trials where no ori change occured
  h2 = plot(-BefSac:AftSac,Auu,'r-');hold on;
  h2b = plot(-BefSac:AftSac,Auu+(2*Asu),'r:');hold on;
  h2b = plot(-BefSac:AftSac,Auu-(2*Asu),'r:');hold on;
  %****** subset of trials where ori change DID occur
  h2 = plot(-BefSac:AftSac,Uuu,'b-');hold on;
  h2b = plot(-BefSac:AftSac,Uuu+(2*Usu),'b:');hold on;
  h2b = plot(-BefSac:AftSac,Uuu-(2*Usu),'b:');hold on;
  %************
  axis tight;
  V = axis;
  axis([-BefSac AftSac 0 V(4)]);
  plot([0,0],[0,V(4)],'k-');
  plot([PreWin(1),PreWin(1)],[0,V(4)],'m-');
  plot([PreWin(2),PreWin(2)],[0,V(4)],'m-');
  xlabel('Time (ms)');
  ylabel('Rate (hz)');
  title('PSTH');
   
  
  %******** we need some plot of saccade end points
  subplot('position',[0.75 0.70 0.19 0.28]);
%   for k = 2:size(TPosList,1)
%       tg = TPosList(k,1);
%       hh = plot(TPosList((k-1):k,2),TPosList((k-1):k,3),[colo(tg),':']);
%       set(hh,'Linewidth',0.5);
%       hold on;
%   end

  colo = 'rkbmgc';
  for k = 1:size(TPosList,1)
     tg = TPosList(k,1);
     h = plot(TPosList(k,2),TPosList(k,3),[colo(tg),'.']); hold on;
     set(h,'Markersize',4);
  end
  plot(0,0,'k+');
  maxo = max(max(abs(TPosList(:,2:3))));
  axis([-maxo maxo -maxo maxo]);
  
  
  %********* plot stim-lock activity 
  ovals = 0:(OriNum-1);
  OriVals = Orivals;
  oset = 1:OriNum;
  
  %***************************
  disp('Running Von Mises Fit');
  JackN = 10;
  NA = size(TrA,1);
  tuneA = [];
  for ik = 1:JackN
      ia = 1 + (ik-1)*floor(NA/JackN);
      ib = ik * floor(NA/JackN);
      subset = [1:ia, ib:NA];
      fitA = vonMises.fit_vonmises(zPreOri(TrA(subset)), StimSpk{1}(TrA(subset)), 0);
      for k = 1:size(OriVals,1)
         tuna(k) = vonMises.vonmises(OriVals(k)*(pi/180),fitA.paramsML);
      end
      tuneA = [tuneA ; tuna];
      disp(sprintf('Computing Jacknife %d of %d',ik,JackN));
  end
  mtuneA = mean(tuneA);
  stuneA = std(tuneA) * sqrt(JackN-1);
  %**********************
  JackN = 10;
  NU = size(TrU,1);
  tuneU = [];
  for ik = 1:JackN
      ia = 1 + (ik-1)*floor(NU/JackN);
      ib = ik * floor(NU/JackN);
      subset = [1:ia, ib:NU];
      fitU = vonMises.fit_vonmises(zPreOri(TrU(subset)), StimSpk{1}(TrU(subset)), 0);
      for k = 1:size(OriVals,1)
         tuna(k) = vonMises.vonmises(OriVals(k)*(pi/180),fitU.paramsML);
      end
      tuneU = [tuneU ; tuna];
      disp(sprintf('Computing Jacknife %d of %d',ik,JackN));
  end
  mtuneU = mean(tuneU);
  stuneU = std(tuneU) * sqrt(JackN-1);
  %************************
  
  subplot('position',[0.75 0.375 0.20 0.25]);
  %******
  for i = 1:OriNum
        zz = find(PreOri == i); % find subset of that orientation
        %****** plot raw counts
        % first plot cases where saccade to RF
        pzz = find( (zPreOri == OriVals(oset(i))) & (TargList == 1) );
        xjit = (rand(size(zPreOri(pzz)))-0.5)*XJit;
        yjit = 1 + (rand(size(StimSpk{1}(pzz)))-0.5)*YJit;      
        H = plot(OriVals(i)*ones(size(pzz))+xjit,StimSpk{1}(pzz) .* yjit,'k.'); hold on;
        set(H,'Markersize',10);
        set(H,'Color',[1,0,0]);
        %****** plot vonMises fit
        h2 = plot(OriVals,mtuneA,'r-');
        set(h2,'Linewidth',2);
        h2 = plot(OriVals,mtuneA + (2 * stuneA),'r-');
        h2 = plot(OriVals,mtuneA - (2 * stuneA),'r-');
        %*********
        % second plot cases where saccade else where
        pzz = find( (zPreOri == OriVals(oset(i))) & (TargList ~= 1) );
        xjit = (rand(size(zPreOri(pzz)))-0.5)*XJit;
        yjit = 1 + (rand(size(StimSpk{1}(pzz)))-0.5)*YJit;      
        H = plot(OriVals(i)*ones(size(pzz))+xjit,StimSpk{1}(pzz) .* yjit,'k.'); hold on;
        set(H,'Markersize',10);
        set(H,'Color',[0,0,1]);
        %****** plot vonMises fit
        h2 = plot(OriVals,mtuneU,'b-');
        set(h2,'Linewidth',2);
        h2 = plot(OriVals,mtuneU + (2 * stuneU),'b-');
        h2 = plot(OriVals,mtuneU - (2 * stuneU),'b-');
        %**************
  end
  axis tight;
  V = axis;
  axis([0 360 0 V(4)]);
  ax = gca;
  set(ax,'XTickLabel',[]);
  ylabel('Spike Counts');
  title('Stim Locked Activity');
  
  %***************************
  disp('Running Von Mises Fit');
  JackN = 10;
  NA = size(TrA,1);
  tuneA = [];
  for ik = 1:JackN
      ia = 1 + (ik-1)*floor(NA/JackN);
      ib = ik * floor(NA/JackN);
      subset = [1:ia, ib:NA];
      fitA = vonMises.fit_vonmises(zPreOri(TrA(subset)), PreSpk{1}(TrA(subset)), 0);
      for k = 1:size(OriVals,1)
         tuna(k) = vonMises.vonmises(OriVals(k)*(pi/180),fitA.paramsML);
      end
      tuneA = [tuneA ; tuna];
      disp(sprintf('Computing Jacknife %d of %d',ik,JackN));
  end
  mtuneA = mean(tuneA);
  stuneA = std(tuneA) * sqrt(JackN-1);
  %**********************
  JackN = 10;
  NU = size(TrU,1);
  tuneU = [];
  for ik = 1:JackN
      ia = 1 + (ik-1)*floor(NU/JackN);
      ib = ik * floor(NU/JackN);
      subset = [1:ia, ib:NU];
      fitU = vonMises.fit_vonmises(zPreOri(TrU(subset)), PreSpk{1}(TrU(subset)), 0);
      for k = 1:size(OriVals,1)
         tuna(k) = vonMises.vonmises(OriVals(k)*(pi/180),fitU.paramsML);
      end
      tuneU = [tuneU ; tuna];
      disp(sprintf('Computing Jacknife %d of %d',ik,JackN));
  end
  mtuneU = mean(tuneU);
  stuneU = std(tuneU) * sqrt(JackN-1);
  %************************
  
  
  subplot('position',[0.75 0.075 0.20 0.25]);
  %*********
  for i = 1:OriNum
        zz = find(PreOri == i); % find subset of that orientation
        %****** plot raw counts
        % first plot cases where saccade to RF
        pzz = find( (zPreOri == OriVals(oset(i))) & (TargList == 1) );
        xjit = (rand(size(zPreOri(pzz)))-0.5)*XJit;
        yjit = 1 + (rand(size(PreSpk{1}(pzz)))-0.5)*YJit;      
        H = plot(OriVals(i)*ones(size(pzz))+xjit,PreSpk{1}(pzz) .* yjit,'k.'); hold on;
        set(H,'Markersize',10);
        set(H,'Color',[1,0,0]);
        %****** plot vonMises fit
        h2 = plot(OriVals,mtuneA,'r-');
        set(h2,'Linewidth',2);
        h2 = plot(OriVals,mtuneA + (2 * stuneA),'r-');
        h2 = plot(OriVals,mtuneA - (2 * stuneA),'r-');
        %*********
        % second plot cases where saccade else where
        pzz = find( (zPreOri == OriVals(oset(i))) & (TargList ~= 1) );
        xjit = (rand(size(zPreOri(pzz)))-0.5)*XJit;
        yjit = 1 + (rand(size(PreSpk{1}(pzz)))-0.5)*YJit;      
        H = plot(OriVals(i)*ones(size(pzz))+xjit,PreSpk{1}(pzz) .* yjit,'k.'); hold on;
        set(H,'Markersize',10);
        set(H,'Color',[0,0,1]);
        %****** plot vonMises fit
        h2 = plot(OriVals,mtuneU,'b-');
        set(h2,'Linewidth',2);
        h2 = plot(OriVals,mtuneU + (2 * stuneU),'b-');
        h2 = plot(OriVals,mtuneU - (2 * stuneU),'b-'); 
        %**************
  end
  axis tight;
  V = axis;
  axis([0 360 0 V(4)]);
  xlabel('Motion direction (degs)');
  ylabel('Spike Counts');
  title('PreSac Locked Activity');
  
  
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
