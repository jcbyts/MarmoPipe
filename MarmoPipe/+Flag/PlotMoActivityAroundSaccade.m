function PlotMoActivityAroundSaccade(Exp,TimeLock,ShowRw,SPClust)
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
  PreSaccade = 1;   % if you only want to view RF responses
  TooFast = 0.10;    % throw out saccades with RT too fast
  %***********************************
  
  if ~isempty(TimeLock)
      TimeLockOff = TimeLock;
  else
      TimeLockOff = 2;  % Time-lock : 0 - saccade onset, 1 - saccade offset
                        %             2 - stim onset
  end
  %******** Parameters for windows on spike counts per trial
  BefSac = 200;  % these are in ms, but below seconds are the spk and eye times
  AftSac = 400;  % so will convert to ms at some point
  %****** selects windows for analysis (non-overlapping now)
  if (1)
     PreWin = [-100,0];
     PostWin = [0,100,200];      
  else
     PreWin = [-300,-200,-100,0];
     PostWin = [0,100,200,300,400,500];
  end
  XJit = 15; %7.5;  % orientation jit
  YJit = 0.1; % percentage rate
  %********************************
  if ~isempty(ShowRw)
    ShowRaw = ShowRw;
  else
    ShowRaw = 0;   % if 1 show raw counts, else show mean curves
  end
  
  %**************
  sp = Exp.sp{SPClust};
  
  %****** experiment specific saccade flagging (if needed)
  SacNum = 0;
  CorNum = 0;
  FixNum = 0;
  Tlist = [];  % list of trials with saccades correct
  TargList = [];  % list of which target he made a saccade towards
  for i = 1:size(Exp.D,1)
        if  strcmp(Exp.D{i}.PR.name,'FlagMo')
           if (Exp.D{i}.PR.error == 0) || (Exp.D{i}.PR.error == 4) % only correct saccade trials (for now) 
               eyeSmo = Exp.D{i}.eyeSmo;
               slist = Exp.D{i}.slist;
               FixRadius = Exp.D{i}.P.fixWinRadius;  % within this bound for be fixating 
               
               %***** this is the cued location (random really, no?)
               % Tx = Exp.D{i}.P.xDeg;
               % Ty = Exp.D{i}.P.yDeg;
               
               %***** instead use the RF location??
               %Tx = Exp.D{i}.P.RF_X;
               %Ty = Exp.D{i}.P.RF_Y;
               
               Tx = Exp.D{i}.PR.targ_x;
               Ty = Exp.D{i}.PR.targ_y;
               TargDist = norm([Tx(1),Ty(1)]);  % eccentricity of target
               BestTarg = NaN;
               
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
                      if (fdist < FixRadius)
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
                               BestTarg = zz(1);
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
                     msactime = NaN
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
                                   end
                                 end
                             end
                         end
                         % disp(sprintf('%d %5.2f %5.2f %5.3f',SacNum,ontime,offstim,(offstim-ontime)));
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
  disp(sprintf('Fixations %d, Correct %d, Saccades Flagged %d',FixNum,CorNum,SacNum)); 
  %*********************************************
  disp(sprintf('Presaccadic, %d att and %d unatt',size(find(TargList == 1),1),...
                                                  size(find(TargList ~= 1),1)));
  
  %******** loop through all trials and build a raster plot time-locked
  %******** on the saccade offset (should see some visual response??)
  Rast = zeros(SacNum,(BefSac+AftSac+1));  % all zeros and ones
  RastList = [];
  SacDur = [];
  MSacDur = [];
  StimOff = [];
  %**** trial by trial stats and counts
  ChangList = [];
  ChangOri = [];
  PreOri = [];
  PostOri = [];
  zPreOri = [];
  zPostOri = [];
  PreSpk = cell(1,size(PreWin,2)-1);
  PostSpk = cell(1,size(PostWin,2)-1);
  %**********
  tr = Tlist(1);
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
      if (TimeLockOff == 2)
         StimOff = [StimOff ; [((1000*(offstim-onstim))) k]];
      else
         StimOff = [StimOff ; [((1000*(offstim-offtime))) k]]; 
      end
      %**********
      if (TimeLockOff == 2)
        SacDur = [SacDur ; [((1000*(ontime-onstim))) k]];
      else    
        SacDur = [SacDur ; [((1000*(ontime-offtime))) k]];
      end
      %****** compute stim onset *****
      z = find( Exp.D{tr}.eyeData(:,5) == 5);  % first moment Matlab detects leaving fix
      if ~isempty(z)
         msacdur = Exp.D{tr}.eyeData(z(1),1) - Exp.D{tr}.eyeData(1,1);
      else
         msacdur = NaN; 
      end
      %****
      if ~isnan(msacdur)
          if (TimeLockOff == 2)
            MSacDur = [MSacDur ; [((1000*(msacdur-onstim))) k]];
          else    
            MSacDur = [MSacDur ; [((1000*(msacdur-offtime))) k]];
          end
      else
          MSacDur = [MSacDur ; [NaN NaN]];
      end
     
      %****************
      if (TimeLockOff == 1)
        staspk = Exp.D{tr}.START_EPHYS + offtime - (BefSac/1000); % covert to secs
        finspk = Exp.D{tr}.START_EPHYS + offtime + (AftSac/1000);
      else
          if (TimeLockOff == 2)
            staspk = Exp.D{tr}.START_EPHYS + onstim - (BefSac/1000); % covert to secs
            finspk = Exp.D{tr}.START_EPHYS + onstim + (AftSac/1000);  
          else
            staspk = Exp.D{tr}.START_EPHYS + ontime - (BefSac/1000); % covert to secs
            finspk = Exp.D{tr}.START_EPHYS + ontime + (AftSac/1000);  
          end
      end
      %**** get spike times from the interval of analysis
      z = find( (sp.st > staspk) & (sp.st < finspk) );
      if ~isempty(z)
         difftime = ceil( 1000 * (sp.st(z) - staspk ) );  % integers in ms
         Rast(k,difftime) = 1;
         RastList = [RastList ; [(difftime-BefSac) ones(size(difftime))*k]];
      end
      %**********
      disp(sprintf('Processing trial %d of %d',k,SacNum));
      %******* this gets more complicated now, as depends on which target
      if (PreSaccade == 1)
         preori = Exp.D{tr}.PR.targori; % before the saccade its targ ori
         apreori = preori;
         postori = Exp.D{tr}.PR.changori;  % change trials
      else
         if ~isnan(tg) 
            ori = Exp.D{tr}.PR.targ_motion(tg);
            dori = ori + Exp.D{tr}.PR.deltaOri;
            if (dori > 360)
                dori = dori - 360;
            end
            if (dori < 0)
                dori = dori + 360;
            end 
            % preori = ori;  % before the saccade its targ ori
            preori = Exp.D{tr}.PR.targori;  % always in RF
            apreori = ori;
            postori = dori;  % change trials 
         else
            preori = NaN;
            apreori = NaN;
            postori = NaN;
         end
      end
      if (Exp.D{tr}.PR.DropStim == 2) % blank trials
         ChangList = [ChangList ; 2];
         ChangOri = [ChangOri ; NaN]; % change to blank
      else
         if (postori ~= apreori)    % orientation change trials
           ChangList = [ChangList ; 1];
           diffori = abs(postori-apreori);
           if (diffori > 180)
            diffori = 360-diffori;
           end
           ChangOri = [ChangOri ; diffori];
         else
           ChangList = [ChangList ; 0];
           ChangOri = [ChangOri ; 0];
         end
      end
      %******** get spike counts pre and post -saccade
      zPreOri =  [zPreOri  ; preori ];  % assumes 12 orientations
      zPostOri = [zPostOri ; postori];  % assumes 12 orientations
      %*** pre-saccadic interval
      for zk = 1:(size(PreWin,2)-1)
          PreSta = PreWin(zk);
          PreEn = PreWin(zk+1);
          onstaspk = Exp.D{tr}.START_EPHYS + ontime + (PreSta/1000); % covert to secs
          onfinspk = Exp.D{tr}.START_EPHYS + ontime + (PreEn/1000);
          z = find( (sp.st >= onstaspk) & (sp.st < onfinspk) );
          if ~isempty(z)
            PreSpk{zk} = [PreSpk{zk} ; (length(z)/(onfinspk-onstaspk)) ];
          else
            PreSpk{zk} = [PreSpk{zk} ; 0];
          end
      end
      %*** post-saccadic interval
      for zk = 1:(size(PostWin,2)-1)
          PostSta = PostWin(zk);
          PostEn = PostWin(zk+1);
          
          %offstaspk = Exp.D{tr}.START_EPHYS + offtime + (PostSta/1000); % covert to secs
          %offfinspk = Exp.D{tr}.START_EPHYS + offtime + (PostEn/1000);
          
          offstaspk = Exp.D{tr}.START_EPHYS + ontime + (PostSta/1000); % covert to secs
          offfinspk = Exp.D{tr}.START_EPHYS + ontime + (PostEn/1000);
          
          z = find( (sp.st >= offstaspk) & (sp.st < offfinspk) );
          if ~isempty(z)
             PostSpk{zk} = [PostSpk{zk} ; (length(z)/(offfinspk-offstaspk)) ];
          else
             PostSpk{zk} = [PostSpk{zk} ; 0];
          end 
      end
      %*******************
  end
  %******** now figure out what is the unique set of orientations
  iz = ~isnan(zPostOri);
  Orivals = unique(zPostOri(iz));
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
  PostOri = [];
  for k = 1:size(zPostOri,1)
      po_ori = find( zPostOri(k) == Orivals ); 
      if ~isempty(po_ori)
         PostOri =  [PostOri  ; po_ori(1)];  % assumes 12 orientations
      else
         PostOri = [PostOri ; NaN];
      end
  end
  %***********
  TrA = find(ChangList == 0);
  TrU = find(ChangList == 1);
  TrB = find(ChangList == 2);
  size(TrA)
  size(TrU)
  size(TrB)
  %****************
  
  %****** plot the raster ****************
  %IH = figure;
  
  H = figure;
  set(H,'Position',[100 100 800 800]);
  subplot('position',[0.1 0.4 0.4 0.5]);
 
  % imagesc(Rast); hold on;
  %******
  % hh = plot(RastList(:,1),RastList(:,2),'k.'); hold on;
  % set(hh,'Markersize',0.2);
  
  %****
  za = find( ChangList == 0);
  zaa = ismember(RastList(:,2),za);
  hh = plot(RastList(zaa,1),RastList(zaa,2),'r.'); hold on;
  set(hh,'Markersize',0.2);
  %********
  zb = find( ChangList == 1);
  zbb = ismember(RastList(:,2),zb);
  hh = plot(RastList(zbb,1),RastList(zbb,2),'b.'); hold on;
  set(hh,'Markersize',0.2);
  %******
  zc = find( ChangList == 2);  % blank trials
  zcc = ismember(RastList(:,2),zc);
  hh = plot(RastList(zcc,1),RastList(zcc,2),'g.'); hold on;
  set(hh,'Markersize',0.2);
  %*******
  
  %**********
  if (TimeLockOff >= 1)
    hh = plot(SacDur(:,1),SacDur(:,2),'k.'); 
    hh2 = plot(MSacDur(:,1),MSacDur(:,2),'b.');
    hh3 = plot(StimOff(:,1),StimOff(:,2),'r.');
  else
    hh = plot(-SacDur(:,1),SacDur(:,2),'k.');      % leave in vpx
    hh2 = plot(-MSacDur(:,1),MSacDur(:,2),'b.');   % leave fix in matlab
    hh3 = plot(-StimOff(:,1),StimOff(:,2),'r.');
  end
  
%   figure(IH);
%   plot(SacDur(:,2),(MSacDur(:,1)-SacDur(:,1)),'k.');
%   xlabel('trial');
%   ylabel('Msac - Vpx sac');
%   figure(H);
  
  %*********
  %set(hh,'Markersize',0.5);
  axis([-BefSac AftSac 0 SacNum]);
  plot([0,0],[0,SacNum],'b-');
  xlabel('Time (ms)');
  ylabel('Trials');
  if (TimeLockOff >= 1)
    if (TimeLockOff == 2)
       title('Saccade Activity (Timelock Stim Onset)'); 
    else
       title('Saccade Activity (Timelock Offset)');
    end
  else
    title('Saccade Activity (Timelock Onset)');   
  end
  %***************************************
  [uu,su] = compute_psth(Rast,Smoothing);
  [Auu,Asu] = compute_psth(Rast(TrA,:),Smoothing);
  [Uuu,Usu] = compute_psth(Rast(TrU,:),Smoothing);
  [Buu,Bsu] = compute_psth(Rast(TrB,:),Smoothing);
  
  %***** plot the psth *******************
  hold on;
  subplot('position',[0.1 0.05 0.4 0.25])
  h2 = plot(-BefSac:AftSac,uu,'k-');hold on;
  set(h2,'Linewidth',2);
  h2b = plot(-BefSac:AftSac,uu+(2*su),'k-');hold on;
  h2b = plot(-BefSac:AftSac,uu-(2*su),'k-');hold on;
  %****** subset of trials where no ori change occured
  h2 = plot(-BefSac:AftSac,Auu,'r-');hold on;
  h2b = plot(-BefSac:AftSac,Auu+(2*Asu),'r:');hold on;
  h2b = plot(-BefSac:AftSac,Auu-(2*Asu),'r:');hold on;
  %****** subset of trials where ori change DID occur
  h2 = plot(-BefSac:AftSac,Uuu,'b-');hold on;
  h2b = plot(-BefSac:AftSac,Uuu+(2*Usu),'b:');hold on;
  h2b = plot(-BefSac:AftSac,Uuu-(2*Usu),'b:');hold on;
  %****** subset of trials where it was blank post sac
  h2 = plot(-BefSac:AftSac,Buu,'g-');hold on;
  h2b = plot(-BefSac:AftSac,Buu+(2*Bsu),'g:');hold on;
  h2b = plot(-BefSac:AftSac,Buu-(2*Bsu),'g:');hold on;
  %************
  axis tight;
  V = axis;
  axis([-BefSac AftSac 0 V(4)]);
  plot([0,0],[0,150],'k-');
  xlabel('Time (ms)');
  ylabel('Rate (hz)');
  title('PSTH for saccade activity');
   
  %*********************
  ovals = 0:(OriNum-1);
  OriVals = Orivals;
  
  %********** Pre-saccadic orientation plot analysis ********** 
  subplot('position',[0.6 0.7 0.3 0.2]);
  % hfz = figure;
  % set(hfz,'position',[400 400 600 500]);
  wsmo = 2;
  if (ShowRaw)
    for zk = 1:(size(PreWin,2)-1)  
      xjit = (rand(size(zPreOri))-0.5)*XJit;
      yjit = 1 + (rand(size(PreSpk{zk}))-0.5)*YJit;
      H = plot(zPreOri+xjit,PreSpk{zk} .* yjit,'ko'); hold on;
      set(H,'Markersize',2);
    end
    axis tight;
  else
    %******** use Jake's fit function instead of other
    if (1)  % dont use, looks to be very slow ... need GPU??
      SampleData = [(TargList==1) zPreOri PreSpk{1}];
      %size(SampleDataC)
      %save('SampleDataC','SampleDataC');
      vonMises.FitPreSaccadicData(SampleData);
    else
            tiny = 1;
            %***********
            for zk = 1:(size(PreWin,2)-1)
              %******
              auu = []; buu = [];
              asulo = []; bsulo = [];
              asuhi = []; bsuhi = [];
              aor = []; bor = [];
              %*******
              % oset = 1:OriNum;
              oset = [floor(OriNum/2):OriNum,1:(floor(OriNum/2)-1)];
              %******
              for i = 1:OriNum
                zz = find(PreOri == i); % find subset of that orientation
                %****** plot raw counts
                % first plot cases where saccade to RF
                pzz = find( (zPreOri == OriVals(oset(i))) & (TargList == 1) );
                xjit = (rand(size(zPreOri(pzz)))-0.5)*XJit;
                yjit = 1 + (rand(size(PreSpk{zk}(pzz)))-0.5)*YJit;      
                H = plot(OriVals(i)*ones(size(pzz))+xjit,PreSpk{zk}(pzz) .* yjit,'k.'); hold on;
                set(H,'Markersize',10);
                set(H,'Color',[1,0,0]);
                %*********
                if ~isempty(pzz)
                   avg = nanmedian(PreSpk{zk}(pzz));
                   savglo = prctile(PreSpk{zk}(pzz),20); %/sqrt(size(PreSpk{zk}(pzz),1));
                   savghi = prctile(PreSpk{zk}(pzz),80);
                   auu = [auu avg];
                   asulo = [asulo savglo];
                   asuhi = [asuhi savghi];
                   aor = [aor OriVals(i)];
                end  
                %*********
                % second plot cases where saccade else where
                pzz = find( (zPreOri == OriVals(oset(i))) & (TargList ~= 1) );
                xjit = (rand(size(zPreOri(pzz)))-0.5)*XJit;
                yjit = 1 + (rand(size(PreSpk{zk}(pzz)))-0.5)*YJit;      
                H = plot(OriVals(i)*ones(size(pzz))+xjit,PreSpk{zk}(pzz) .* yjit,'k.'); hold on;
                set(H,'Markersize',10);
                set(H,'Color',[0,0,1]);
                %**************
                if ~isempty(pzz)
                   avg = nanmedian(PreSpk{zk}(pzz));
                   savglo = prctile(PreSpk{zk}(pzz),20); %/sqrt(size(PreSpk{zk}(pzz),1));
                   savghi = prctile(PreSpk{zk}(pzz),80);
                   buu = [buu avg];
                   bsulo = [bsulo savglo];
                   bsuhi = [bsuhi savghi];
                   bor = [bor OriVals(i)];
                end
              end
              %*****
              wsmo = 2;
              auu = tunesmooth(auu,wsmo);
              asulo = tunesmooth(asulo,wsmo);
              asuhi = tunesmooth(asuhi,wsmo);
              H = plot(aor,auu,'r-'); hold on;
              set(H,'LineWidth',2);
              H = plot(aor,asulo,'r-'); hold on;
              H = plot(aor,asuhi,'r-'); hold on;
              %*****
              buu = tunesmooth(buu,wsmo);
              bsulo = tunesmooth(bsulo,wsmo);
              bsuhi = tunesmooth(bsuhi,wsmo);
              H = plot(bor,buu,'b-'); hold on;
              set(H,'LineWidth',2);
              H = plot(bor,bsulo,'b-'); hold on;
              H = plot(bor,bsuhi,'b-'); hold on;
              % ******
            end
            axis tight;
            V = axis;
            axis([V(1) V(2) 0 V(4)]);
            endx = V(1) + 0.7*(V(2)-V(1));
            endy = 0.85*V(4);
            h2 = text(endx,endy,'Towards RF');
            set(h2,'Color',[1,0,0]);
            endx = V(1) + 0.7*(V(2)-V(1));
            endy = 0.75*V(4);
            h2 = text(endx,endy,' Away RF');
            set(h2,'Color',[0,0,1]);

          end
          xlabel('Direction (deg)');
          ylabel('Firing (sp/s)');
          title('Pre-saccadic Tuning');
  end
  
  %********* Analyze Post-saccadic where orientation does not change***
  subplot('position',[0.6 0.4 0.3 0.2]);
  if (ShowRaw)
    z = find( ChangList == 0);
    for zk = 1:(size(PostWin,2)-1)
       xjit = (rand(size(zPostOri(z)))-0.5)*XJit;
       yjit = 1 + (rand(size(PostSpk{zk}(z)))-0.5)*YJit;      
       H = plot(zPostOri(z)+xjit,PostSpk{zk}(z) .* yjit,'ro'); hold on;
       set(H,'Markersize',2);
    end
    axis tight;
  else
    %*********
    iAuu = [];
    for zk = 1:(size(PostWin,2)-1)
        z = find( ChangList == 0);
        %************
        Auu = [];
        Asu = [];
        Aor = [];
        for i = 1:OriNum
           zz = find( (PostOri== i) & (ChangList == 0) ); % find subset of that orientation
           %****** plot raw counts
           pzz = find(zPostOri(zz) == OriVals(i));
           xjit = (rand(size(zPostOri(zz(pzz))))-0.5)*XJit;
           yjit = 1 + (rand(size(PostSpk{zk}(zz(pzz))))-0.5)*YJit;      
           H = plot(zPostOri(zz(pzz))+xjit,PostSpk{zk}(zz(pzz)) .* yjit,'ro'); hold on;
           set(H,'Markersize',2);
           set(H,'Color',[1,0.5,0.5]); 
           %**************
           if ~isempty(zz)
             avg = nanmean(PostSpk{zk}(zz));
             savg = nanstd(PostSpk{zk}(zz))/sqrt(size(PostSpk{zk}(zz),1));
             Auu = [Auu avg];
             Asu = [Asu savg];
             Aor = [Aor OriVals(i)];
           end
        end
        iAuu = [iAuu ; Auu];
    end
    Auu = mean(iAuu);
    Asu = std(iAuu)/sqrt(size(iAuu,1));
    H = plot(Aor,Auu,'r.-');
    set(H,'Linewidth',2);
    plot(Aor,(Auu+(2*Asu)),'r-');
    plot(Aor,(Auu-(2*Asu)),'r-');    
    axis tight;
  end  
  %********** Analyze Post-saccadic where orientation does change
  if (ShowRaw)
    z = find( ChangList == 1);
    for zk = 1:(size(PostWin,2)-1)
        xjit = (rand(size(zPostOri(z)))-0.5)*XJit;
        yjit = 1 + (rand(size(PostSpk{zk}(z)))-0.5)*YJit;      
        H = plot(zPostOri(z)+xjit,PostSpk{zk}(z) .* yjit,'bo'); hold on;
        set(H,'Markersize',2);
    end
    axis tight;
  else
    %**********
    z = find( ChangList == 1);
    %**********
    iUuu = [];
    for zk = 1:(size(PostWin,2)-1)
        %**********
        Uuu = [];
        Usu = [];
        Uor = [];
        for i = 1:OriNum
           zz = find( (PostOri == i) & (ChangList == 1) ); % find subset of that orientation
           %****** plot raw counts
           pzz = find(zPostOri(zz) == OriVals(i));
           xjit = (rand(size(zPostOri(zz(pzz))))-0.5)*XJit;
           yjit = 1 + (rand(size(PostSpk{zk}(zz(pzz))))-0.5)*YJit;      
           H = plot(zPostOri(zz(pzz))+xjit,PostSpk{zk}(zz(pzz)) .* yjit,'bo'); hold on;
           set(H,'Markersize',2);
           set(H,'Color',[0.5,0.5,1]); 
           %**************       
           if ~isempty(zz)
               avg = nanmean(PostSpk{zk}(zz));
               savg = nanstd(PostSpk{zk}(zz))/sqrt(size(PostSpk{zk}(zz),1));
               Uuu = [Uuu avg];
               Usu = [Usu savg];
               Uor = [Uor OriVals(i)];
           end
        end
        iUuu = [iUuu ; Uuu];
    end
    Uuu = mean(iUuu);
    Usu = std(iUuu)/sqrt(size(iUuu,1));
    H = plot(Uor,Uuu,'b.-');
    set(H,'Linewidth',2);
    plot(Uor,(Uuu+(2*Usu)),'b-');
    plot(Uor,(Uuu-(2*Usu)),'b-');
    axis tight;    
  end
  xlabel('Orientation');
  ylabel('Firing (sp/s)');
  title('Post-saccadic Ori (Change)');
  %******************************
  
  subplot('position',[0.6 0.1 0.3 0.2]);
  if (ShowRaw)
    z = find( ChangList == 2);
    for zk = 1:(size(PostWin,2)-1)
      xjit = (rand(size(zPreOri(z)))-0.5)*XJit;
      yjit = 1 + (rand(size(PostSpk{zk}(z)))-0.5)*YJit;      
      H = plot(zPreOri(z)+xjit,PostSpk{zk}(z) .* yjit,'go'); hold on;
      set(H,'Markersize',2);
    end
    axis tight;
  else
    %*********  
    z = find( ChangList == 2);
    iBuu = [];
    for zk = 1:(size(PostWin,2)-1)
        %***********
        Buu = [];
        Bsu = [];
        Bor = [];
        for i = 1:OriNum
           zz = find( (PreOri== i) & (ChangList == 2) ); % find subset of that orientation
           %****** plot raw counts
           pzz = find(zPreOri(zz) == OriVals(i));
           xjit = (rand(size(zPreOri(zz(pzz))))-0.5)*XJit;
           yjit = 1 + (rand(size(PostSpk{zk}(zz(pzz))))-0.5)*YJit;      
           H = plot(zPreOri(zz(pzz))+xjit,PostSpk{zk}(zz(pzz)) .* yjit,'go'); hold on;
           set(H,'Markersize',2);
           set(H,'Color',[0.5,1,0.5]); 
           %**************       
           if ~isempty(zz)
             avg = nanmean(PostSpk{zk}(zz));
             savg = nanstd(PostSpk{zk}(zz))/sqrt(size(PostSpk{zk}(zz),1));
             Buu = [Buu avg];
             Bsu = [Bsu savg];
             Bor = [Bor OriVals(i)];
           end
        end
        iBuu = [iBuu ; Buu];
    end
    Buu = mean(iBuu);
    Bsu = std(iBuu)/sqrt(size(iBuu,1));
    H = plot(Bor,Buu,'g.-');
    set(H,'Linewidth',2);
    plot(Bor,(Buu+(2*Bsu)),'g-');
    plot(Bor,(Buu-(2*Bsu)),'g-');
    axis tight;  
  end  
  
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

function y = tunesmooth(x,sig)
   N = length(x);
   y = x;
   for k = 1:N
      sumo = 0;
      wsumo = 0;
      for j = 1:N
         dist = abs(j-k);
         if (dist > (N/2))
             dist = N - dist;
         end
         val = exp(-0.5*(dist/sig)^2);
         sumo = sumo + (val * x(j));
         wsumo = wsumo + val;
      end
      if (wsumo)
        y(k) = (sumo/wsumo);
      end
   end
return;