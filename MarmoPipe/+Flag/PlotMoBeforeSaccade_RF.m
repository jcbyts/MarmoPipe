function Info = PlotMoActivityAroundSaccade_RF(Exp,SPClust,StimX,StimY,GRID,FileTag)
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
%     Info - a struct with necessary results for plots and for pop analyses

  %*** it will go through all Camoflag trials and determine 
  %*** the saccades into the aperture
  EyeRadius = 2.5;  % should be defined in protocol, but what we count
                    % as the saccade end point landing close enough to 
                    % the target
  AmpThresh = 0.4;  % fraction of saccade distance of total target ecc
                    % before we include that saccade as totally correct
  Smoothing = 10.0; % Gauss smoothing for PSTH plot
  TooFast = 0.10;   % throw out saccades with RT too fast
  JackN = 10;       % Jacknife 
  %***********************************
  
  %******** Parameters for windows on spike counts per trial
  BefSac = 300;  % these are in ms, but below seconds are the spk and eye times
  AftSac = 50;  % so will convert to ms at some point
  BefStim = 50;  % these are in ms, but below seconds are the spk and eye times
  AftStim = 300;  % so will convert to ms at some point
  %****** selects windows for analysis (non-overlapping now)
  PreWin = [-100,0]; %[-50,25];%[-100,0];% for plotting purposes
  StimWin = [50,150]; %[50,150];%[50,170]; % for plotting
  %*********
  PreWinWid = 40;
  PreWinBin = [-BefSac:10:(AftSac-PreWinWid)];
  NPreBin = size(PreWinBin,2)+1;
  StimWinWid = 40;
  StimWinBin = [-BefStim:10:(AftStim-StimWinWid)];
  NStimBin = size(StimWinBin,2)+1;  
  %**********
  XJit = 10; % orientation jit
  YJit = 0.05; % percentage rate
  %********************************
  
  %****** store Info here ***********
  Info.FileTag = FileTag;
  Info.SPClust = SPClust;
  Info.GRID = GRID;
  Info.EyeRadius = EyeRadius;
  Info.AmpThresh = AmpThresh;
  Info.Smoothing = Smoothing;
  Info.TooFast = TooFast;   
  Info.JackN = JackN;
  Info.BefSac = BefSac;
  Info.AftSac = AftSac;  
  Info.BefStim = BefStim; 
  Info.AftStim = AftStim; 
  Info.PreWin = PreWin;
  Info.StimWin = StimWin;
  Info.PreWinWid = PreWinWid;
  Info.PreWinBin = PreWinBin;
  Info.NPreBin = NPreBin;
  Info.StimWinWid = StimWinWid;
  Info.StimWinBin = StimWinBin;
  Info.NStimBin = NStimBin;
  Info.XJit = XJit;
  Info.YJit = YJit;
  %**********************************
  
  %**************
  sp = Exp.sp{SPClust};
  
  %****** experiment specific saccade flagging (if needed)
  SacNum = 0;
  CorNum = 0;
  FixNum = 0;
  Tlist = [];  % list of trials with saccades correct
  TargList = [];  % list of which target he made a saccade towards
  TPosList = [];  % list of target choice, and x,y of sac endpoint
  TPosTrace = cell(1,1);
  PosCnt = 1;
  SecondSac = 0;
  for i = 1:size(Exp.D,1)
        if  strcmp(Exp.D{i}.PR.name(1:6),'FlagMo')
         
            %if (1) % (Exp.D{i}.PR.deltaOri == 0)  % isnan(Exp.D{i}.PR.changori)  
         
         if isfield(Exp.D{i}.PR,'catchtrial')
             if (Exp.D{i}.PR.catchtrial == 1)  % do not include this trial
                 continue;
             end
         end
         
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
                  %***** look for any 2nd saccade before the big one
                  
                  %****************
                  if ~isempty(offsac)
                      ontime = slist(offsac(1),1);  % onset of that saccade
                      if ( (offstim-ontime) < 0.03)
                         if ~isnan(msactime) && (msactime > ontime)     
                             FixNum = FixNum + foundfix;
                             CorNum = CorNum + foundcor;         
                             SacNum = SacNum + foundsac;
                             %****** final check, is there a 2nd saccade
                             sac2nd = find( (slist(:,2) >= onstim) & (slist(:,2) < eyeSmo(stasac,1)) );
                             if (~isempty(sac2nd))
                                 SecondSac = SecondSac + 1;
                                 foundsac = 0;
                             end
                             %********
                             if (foundsac == 1) 
                                 if ~isnan(msactime) && ~isnan(onstim)
                                   if ( (msactime-onstim) > TooFast)
                                     if isfield(Exp.D{i},'START_EPHYS')  
                                       Tlist = [Tlist ; i];
                                       TargList = [TargList ; BestTarg];
                                       %*********************
                                       TPosList = [TPosList ; [BestTarg,BestEtx,BestEty]];
                                       %************
                                       tzz = find( (eyeSmo(:,1) >= onstim) & ...
                                                  (eyeSmo(:,1) < (eyeSmo(stasac,1)) ) );
                                       TPosTrace{PosCnt} = eyeSmo(tzz,1:3);
                                       PosCnt = PosCnt + 1;
                                       %************
                                     end
                                   end
                                 end
                             end
                             %********
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
  %*************************
  disp(' ');
  disp(sprintf('Fixations %d, Correct %d', FixNum,CorNum));
  disp(sprintf('Saccades Flagged %d SecondSac %d',SacNum,SecondSac));
  disp('  ');
  %*********************************************
  disp(sprintf('Presaccadic, %d att and %d unatt',size(find(TargList == 1),1),...
                                                  size(find(TargList ~= 1),1)));
                                              
  %**********
  Info.FixNum = FixNum;
  Info.CorNum = CorNum;
  Info.SacNum = SacNum;
  Info.SecondSac = SecondSac;
  Info.ANum = size(find(TargList == 1));
  Info.UNum = size(find(TargList ~= 1));
  %****************
  
  
  %******** loop through all trials and build a raster plot time-locked
  %******** on the saccade offset (should see some visual response??)
  StimRast = zeros(SacNum,(BefStim+AftStim+1));  % all zeros and ones
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
  PreSpk = cell(1,NPreBin);
  StimSpk = cell(1,NStimBin); 
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
            staspk = Exp.D{tr}.START_EPHYS + onstim - (BefStim/1000); % covert to secs
            finspk = Exp.D{tr}.START_EPHYS + onstim + (AftStim/1000);  
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
                 StimRastList = [StimRastList ; [(difftime-BefStim) ones(size(difftime))*k]];
             else
                 SacRast(k,difftime) = 1;
                 SacRastList = [SacRastList ; [(difftime-BefSac) ones(size(difftime))*k]];                    
             end
         end

         if (TimeLockOff == 2)
            for zk = 1:NStimBin 
              if (zk == NStimBin)
                PreSta = StimWin(1);
                PreEn = StimWin(2);
              else
                PreSta = StimWinBin(zk);
                PreEn = StimWinBin(zk) + StimWinWid;     
              end
              onstaspk = Exp.D{tr}.START_EPHYS + onstim + (PreSta/1000); % covert to secs
              onfinspk = Exp.D{tr}.START_EPHYS + onstim + (PreEn/1000);
              z = find( (sp.st >= onstaspk) & (sp.st < onfinspk) );
              if ~isempty(z)
                 StimSpk{zk} = [StimSpk{zk} ; (length(z)/(onfinspk-onstaspk)) ];
              else
                 StimSpk{zk} = [StimSpk{zk} ; 0];
              end
            end
         else
           for zk = 1:NPreBin
              if (zk == NPreBin) 
                PreSta = PreWin(1);
                PreEn = PreWin(2);
              else
                PreSta = PreWinBin(zk);
                PreEn = PreSta + PreWinWid;
              end
              onstaspk = Exp.D{tr}.START_EPHYS + ontime + (PreSta/1000); % covert to secs
              onfinspk = Exp.D{tr}.START_EPHYS + ontime + (PreEn/1000);
              z = find( (sp.st >= onstaspk) & (sp.st < onfinspk) );
              if ~isempty(z)
                 PreSpk{zk} = [PreSpk{zk} ; (length(z)/(onfinspk-onstaspk)) ];
              else
                 PreSpk{zk} = [PreSpk{zk} ; 0];
              end
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
  OriVals = Orivals;
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
  TrU = find(ChangList ~= 1);
  disp(sprintf('Into RF %d and Away %d',size(TrA,1),size(TrU,1)));
  %****************
  
  %************
  Info.Orivals = Orivals;
  Info.OriVals = Orivals;
  Info.OriNum = OriNum;
  Info.TrA = TrA;
  Info.TrU = TrU;
  Info.PreOri = PreOri;
  Info.zPreOri = zPreOri;
  Info.PreSpk = PreSpk;
  Info.StimSpk = StimSpk;
  Info.tr = tr;
  Info.ChangList = ChangList;
  %*****
  Info.StimRast = StimRast;
  Info.StimRastList = StimRastList;
  Info.SacRast = SacRast;
  Info.SacRastList = SacRastList;
  %***********
  Info.StimSacDur = StimSacDur;
  Info.StimMSacDur = StimMSacDur;
  Info.SacSacDur = SacSacDur;
  Info.SacMSacDur = SacMSacDur;
  %*****************
  
  %****** plot the Spatial RF here ***********
  if ~isempty(StimX) && ~isempty(StimY)
     [RFsets,tXX,Zx,Zy,RFmino,RFmaxo] = Forage.PlotForageSpatialKernel_RF(StimX,StimY,GRID); 
      %********
      Info.RFsets = RFsets;
      Info.tXX = tXX;
      Info.Zx = Zx;
      Info.Zy = Zy;
      Info.RFmino = RFmino;
      Info.RFmaxo = RFmaxo;
  else
      Info.RFsets = [];
      Info.tXX = [];
      Info.Zx = [];
      Info.Zy = [];
      Info.RFmino = [];
      Info.RFmaxo = [];    
  end
  %*******************
  
  %************* compute PSTHs for stim locked
  [AuuStim,AsuStim] = compute_psth(StimRast(TrA,:),zPreOri(TrA),OriVals,Smoothing);
  [UuuStim,UsuStim] = compute_psth(StimRast(TrU,:),zPreOri(TrU),OriVals,Smoothing); 
  %***************************************
  [AuuPre,AsuPre] = compute_psth(SacRast(TrA,:),zPreOri(TrA),OriVals,Smoothing);
  [UuuPre,UsuPre] = compute_psth(SacRast(TrU,:),zPreOri(TrU),OriVals,Smoothing);
  %**********
  Info.AuuStim = AuuStim;
  Info.AsuStim = AsuStim;
  Info.UuuStim = UuuStim;
  Info.UsuStim = UsuStim;
  Info.AuuPre = AuuPre;
  Info.AsuPre = AsuPre;
  Info.UuuPre = UuuPre;
  Info.UsuPre = UsuPre;
  %********************
  
  %******** eye movement analysis of fixation positions
  TN = max(TPosList(:,1));
  TPosTraceList = cell(1,TN);
  for k = 1:size(TPosList,1)
     tg = TPosList(k,1);
     TPosTraceList{1,tg} = [TPosTraceList{1,tg} ; [mean(TPosTrace{k}(:,2)),mean(TPosTrace{k}(:,3))] ];
  end
  %****** get mean eye position before saccade
  TuPos = zeros(TN,4);
  for k = 1:size(TPosTraceList,2)
      TuPos(k,1) = mean(TPosTraceList{1,k}(:,1));
      TuPos(k,2) = mean(TPosTraceList{1,k}(:,2));
      zz = find(TPosList(:,1) == k);
      TuPos(k,3) = mean(TPosList(zz,2));
      TuPos(k,4) = mean(TPosList(zz,3));
  end
  %**************
  Info.TN = TN;
  Info.TPosList = TPosList;
  Info.TPosTraceList = TPosTraceList;
  Info.TPosTrace = TPosTrace;
  Info.TuPos = TuPos;
  %******************

  %******** Compute tuning curve information
  OriVals = Orivals;
  1
  [Afitstim,Atunestim,Awidstim,Arawstim,Afvalstim] = FitWithVonMises(OriVals,JackN,TrA,zPreOri,StimSpk{NStimBin});
  2
  [Ufitstim,Utunestim,Uwidstim,Urawstim,Ufvalstim] = FitWithVonMises(OriVals,JackN,TrU,zPreOri,StimSpk{NStimBin});
  TrT = [TrA ; TrU]; % everything
  3
  [Tfitstim,Ttunestim,Twidstim,Trawstim,Tfvalstim] = FitWithVonMises(OriVals,JackN,TrT,zPreOri,StimSpk{NStimBin});
  %*******
  4
  [Afitpre,Atunepre,Awidpre,Arawpre,Afvalpre] = FitWithVonMises(OriVals,JackN,TrA,zPreOri,PreSpk{NPreBin});
  5
  [Ufitpre,Utunepre,Uwidpre,Urawpre,Ufvalpre] = FitWithVonMises(OriVals,JackN,TrU,zPreOri,PreSpk{NPreBin});
  %**********
  Info.Afitstim = Afitstim;
  Info.Atunestim = Atunestim;
  Info.Awidstim = Awidstim;
  Info.Ufitstim = Ufitstim;
  Info.Utunestim = Utunestim;
  Info.Uwidstim = Uwidstim;
  Info.Afitpre = Afitpre;
  Info.Atunepre = Atunepre;
  Info.Awidpre = Awidpre;
  Info.Ufitpre = Ufitpre;
  Info.Utunepre = Utunepre;
  Info.Uwidpre = Uwidpre;
  Info.TrT = TrT;
  Info.Tfitstim = Tfitstim;
  Info.Ttunestim = Ttunestim;
  Info.Twidstim = Twidstim;
  Info.Arawstim = Arawstim;
  Info.Urawstim = Urawstim;
  Info.Trawstim = Trawstim;
  Info.Arawpre = Arawpre;
  Info.Urawpre = Urawpre;
  Info.Afvalstim = Afvalstim;
  Info.Ufvalstim = Ufvalstim;
  Info.Tfvalstim = Tfvalstim;
  Info.Afvalpre = Afvalpre;
  Info.Ufvalpre = Ufvalpre;
  %*************  
  
  %******* Compute Towards (Att) DSI 
  wsum = 0;
  Wvec = 0;
  for k = 1:size(OriVals,1)
     ori = OriVals(k)*(pi/180);
     vc = complex(cos(ori),sin(ori));
     wsum = wsum + Atunepre.mu(k);
     Wvec = Wvec + (Atunepre.mu(k) * vc);
  end
  if (wsum > 0)
      ADSI = abs( Wvec / wsum );
      Wvec = Wvec / wsum;
  else
      ADSI = 0;
  end
  %******* Compute Towards (Att) DSI 
  wsum = 0;
  Wvec = 0;
  for k = 1:size(OriVals,1)
     ori = OriVals(k)*(pi/180);
     vc = complex(cos(ori),sin(ori));
     wsum = wsum + Utunepre.mu(k);
     Wvec = Wvec + (Utunepre.mu(k) * vc);
  end
  if (wsum > 0)
      UDSI = abs( Wvec / wsum );
      Wvec = Wvec / wsum;
  else
      UDSI = 0;
  end
  %******* Compute DSI and preferred direction
  wsum = 0;
  Wvec = 0;
  for k = 1:size(OriVals,1)
     ori = OriVals(k)*(pi/180);
     vc = complex(cos(ori),sin(ori));
     wsum = wsum + Ttunestim.mu(k);
     Wvec = Wvec + (Ttunestim.mu(k) * vc);
  end
  if (wsum > 0)
      DSI = abs( Wvec / wsum );
      Wvec = Wvec / wsum;
  else
      DSI = 0;
  end
  %******** find preferred direction
  PrefDir = 0;
  dist = -Inf;
  for k = 1:size(OriVals,1)
      ori = OriVals(k)*(pi/180);
      vc = complex(cos(ori),sin(ori));
      dotprod = (real(vc) * real(Wvec)) + (imag(vc) * imag(Wvec));
      if (dotprod > dist)
         PrefDir = k;
         dist = dotprod;
      end
  end
  %***************
  Info.Wvec = Wvec;
  Info.DSI = DSI;
  Info.ADSI = ADSI;
  Info.UDSI = UDSI;
  Info.PrefDir = PrefDir;
  %*******
  NPrefDir = PrefDir + (size(OriVals,1)/2);
  if (NPrefDir > size(OriVals,1))
      NPrefDir = NPrefDir - size(OriVals,1);
  end
  Info.NPrefDir = NPrefDir;
  %*********************
  
  %************* compute PSTHs for stim locked
  PrefWid = 1;
  PrefSet = [];
  NPrefSet = [];
  for k = 1:size(OriVals,1)
      dist = abs(k-PrefDir);
      if (dist > (size(OriVals,1)/2))
         dist = size(OriVals,1) - dist; 
      end
      if (dist <= PrefWid)
          PrefSet = [PrefSet ; k];
      end
  end
  for k = 1:size(OriVals,1)
      dist = abs(k-NPrefDir);
      if (dist > (size(OriVals,1)/2))
         dist = size(OriVals,1) - dist; 
      end
      if (dist <= PrefWid)
          NPrefSet = [NPrefSet ; k];
      end
  end
  %**************
  [AuuStimPref,AsuStimPref] = compute_psth(StimRast(TrA,:),zPreOri(TrA),OriVals(PrefSet),Smoothing);
  [AuuStimNPref,AsuStimNPref] = compute_psth(StimRast(TrA,:),zPreOri(TrA),OriVals(NPrefSet),Smoothing);
  [UuuStimPref,UsuStimPref] = compute_psth(StimRast(TrU,:),zPreOri(TrU),OriVals(PrefSet),Smoothing); 
  [UuuStimNPref,UsuStimNPref] = compute_psth(StimRast(TrU,:),zPreOri(TrU),OriVals(NPrefSet),Smoothing); 
  %***************************************
  [AuuPrePref,AsuPrePref] = compute_psth(SacRast(TrA,:),zPreOri(TrA),OriVals(PrefSet),Smoothing);
  [AuuPreNPref,AsuPreNPref] = compute_psth(SacRast(TrA,:),zPreOri(TrA),OriVals(NPrefSet),Smoothing);
  [UuuPrePref,UsuPrePref] = compute_psth(SacRast(TrU,:),zPreOri(TrU),OriVals(PrefSet),Smoothing);
  [UuuPreNPref,UsuPreNPref] = compute_psth(SacRast(TrU,:),zPreOri(TrU),OriVals(NPrefSet),Smoothing);
  %**********
  Info.PrefSet = PrefSet;
  Info.NPrefSet = NPrefSet;
  %*****
  Info.AuuStimPref = AuuStimPref;
  Info.AsuStimPref = AsuStimPref;
  Info.UuuStimPref = UuuStimPref;
  Info.UsuStimPref = UsuStimPref;
  Info.AuuPrePref = AuuPrePref;
  Info.AsuPrePref = AsuPrePref;
  Info.UuuPrePref = UuuPrePref;
  Info.UsuPrePref = UsuPrePref;
  %*****
  Info.AuuStimNPref = AuuStimNPref;
  Info.AsuStimNPref = AsuStimNPref;
  Info.UuuStimNPref = UuuStimNPref;
  Info.UsuStimNPref = UsuStimNPref;
  Info.AuuPreNPref = AuuPreNPref;
  Info.AsuPreNPref = AsuPreNPref;
  Info.UuuPreNPref = UuuPreNPref;
  Info.UsuPreNPref = UsuPreNPref;
  %********************
  
  % pr_ori = find( zPreOri(k) == Orivals ); 
  
  %******** ROC analysis across the different time bins *********
  ARoc = zeros(NPreBin,3);  %AUC, error bounds: low and high
  URoc = zeros(NPreBin,3);
  TimeRoc = PreWinBin + (PreWinWid/2);
  for zk = 1:NPreBin
      %****
      attpref = [];   % list of counts
      for k = 1:size(PrefSet,1)
          zz = find( zPreOri(TrA) == OriVals(PrefSet(k)) );
          attpref = [attpref ; PreSpk{zk}(TrA(zz))];         
      end
      %********
      attnpref = [];  % list of counts
      for k = 1:size(NPrefSet,1)
          zz = find( zPreOri(TrA) == OriVals(NPrefSet(k)) );
          attnpref = [attnpref ; PreSpk{zk}(TrA(zz))];         
      end
      %********
      auc = Flag.JakeROC([attpref ; attnpref]',...
              [ones(size(attpref)) ; zeros(size(attnpref))]',0.05);
      ARoc(zk,1) = auc.AUC;
      ARoc(zk,2) = auc.ci(1);
      ARoc(zk,3) = auc.ci(2);
      %****************************
      unattpref = [];   % list of counts
      for k = 1:size(PrefSet,1)
          zz = find( zPreOri(TrU) == OriVals(PrefSet(k)) );
          unattpref = [unattpref ; PreSpk{zk}(TrU(zz))];         
      end
      %********
      unattnpref = [];  % list of counts
      for k = 1:size(NPrefSet,1)
          zz = find( zPreOri(TrU) == OriVals(NPrefSet(k)) );
          unattnpref = [unattnpref ; PreSpk{zk}(TrU(zz))];         
      end
      %********
      auc = Flag.JakeROC([unattpref ; unattnpref]',...
              [ones(size(unattpref)) ; zeros(size(unattnpref))]',0.05);
      URoc(zk,1) = auc.AUC;
      URoc(zk,2) = auc.ci(1);
      URoc(zk,3) = auc.ci(2);
      %****************************
  end
  %*********
  Info.ARoc = ARoc;
  Info.URoc = URoc;
  Info.TimeRoc = TimeRoc;
  %**************************************************************
  
  Flag.PlotInfo_MoBeforeSaccade(Info);  % plot out everything
   
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

function [u_smooth,sem_smooth] = compute_psth(Raster,TrialOris,OriVals,Smoothing)
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
    
    traces = [];
    for k = 1:size(OriVals,1)
        zz = find( TrialOris(excludeset) == OriVals(k) );
        if ~isempty(zz)
             if (size(zz,1) == 1)
                 psth = Raster(excludeset(zz),:);
             else
                 psth = mean(Raster(excludeset(zz),:));
             end
             psth = psth*1000;
             traces = [traces ; psth];
        end
    end
    smooth_data2 = gauss_smooth(mean(traces),Smoothing);

    smoothsub = [smoothsub; smooth_data2];
   end
   u_smooth = mean(smoothsub);
   sem_smooth = (JNum-1) * var(smoothsub);
   sem_smooth = sqrt(sem_smooth);   % variance is multiplied by N-1 Jacknife

return;

function [Afit,Atune,Awid,Araw,Fval] = FitWithVonMises(OriVals,JackN,TrA,zPreOri,StimSpk)
  
  %***************************
  disp('Running Von Mises Fit');
  NA = size(TrA,1);
  mfitA = [];
  tuneA = [];
  wfitA = [];
  fvalA = [];
  
  StimSpk
  
  for ik = 1:JackN
      ia = 1 + (ik-1)*floor(NA/JackN);
      ib = ik * floor(NA/JackN);
      subset = [1:ia, ib:NA];
      fitA = vonMises.fit_vonmises(zPreOri(TrA(subset)), StimSpk(TrA(subset)), 0);
      for k = 1:size(OriVals,1)
         tuna(k) = vonMises.vonmises(OriVals(k)*(pi/180),fitA.paramsML);
      end
      tuneA = [tuneA ; tuna];
      mfitA = [mfitA ; fitA.paramsML];
      fvalA = [fvalA ; fitA.fvalue];
      %******* compute half-width
      for k = 1:360
          ztuna(k) = vonMises.vonmises((k*(pi/180)),fitA.paramsML);
      end
      maxo = max(ztuna);
      mino = min(ztuna);
      half = 0.5*(maxo+mino);
      zz = find( ztuna > half);
      wid = length(zz);
      %**************************
      wfitA = [wfitA ; wid];
      disp(sprintf('Computing Jacknife %d of %d',ik,JackN));
  end
  %********
  mtuneA = mean(tuneA);
  stuneA = std(tuneA) * sqrt(JackN-1);
  Atune.mu = mtuneA;
  Atune.sem = stuneA;
  %******
  mfitmu = mean(mfitA);
  mfitsem = std(mfitA) * sqrt(JackN-1);
  Afit.mu = mfitmu; 
  Afit.sem = mfitsem; 
  %*********
  Fval = mean(fvalA);  % return goodness of fit
  %*******
  Awid.mu = mean(wfitA);
  Awid.sem = std(wfitA) * sqrt(JackN-1);
  %******* compute the raw tuning curve as well
  Araw = [];
  for k = 1:size(OriVals,1)
      zz = find( zPreOri(TrA) == OriVals(k) );
      Araw(k) = mean( StimSpk(TrA(zz)));
  end
  %************************
  
 return;
  
