function [SacX] = SacMatrix_ForageGratingKernel(Exp)
% function StimMatrix_ForageGratingKernel(Exp)
%   input: takes Exp struct as the input
%          SPClust - number of spike cluster to analyze
%  Exp fields:
%    D: {100x1 cell}  - data per trial
%    S: [1x1 struct]  - rig settings struct
%   sp: [1x1 struct]  - spikes from spike sorting
%   cp: (will be a 60hz filtered LFP at 1000 hz)
%   ey: (will be a struct for eye position)
% Outputs:
%   SacX:  [Nframes,3] 
%        1st column    - 0-no saccade, 
%                        1- frame of saccade onset to pref target
%                        2- frame of saccade onset to non-pref target
%                        3- other saccade
%        2nd column    - x amp of saccade
%        3rd column    - y amp of saccade
%        4th column    - 0 - no saccade
%                      - 1,2,3 as below, but offset
%*************

    TargWin = 2.0;  % what counts as on the target
    MinSacDur = 0.20; % min saccade duration 
    %**** go through trials, building a stimulus representations
    %**** in the coordinate frame of spike times ***************
    %**** at this time, no consideration of saccades or probe
    %**** stimulus timing is considered (later)
    %***
    %***** first scroll trials to compute size of stim array in time
    Nstim = 0;
    for zk = 1:size(Exp.D,1)
       if (length(Exp.D{zk}.PR.name) >= 6)
        if strcmp(Exp.D{zk}.PR.name(1:6),'Forage')
          if (Exp.D{zk}.PR.noisetype == 1) % Hartley
            Nstim = Nstim + size(Exp.D{zk}.PR.NoiseHistory,1);
          end
        end
       end
    end
    
    SacOnsets = [];
    SacOffsets = [];
    StimOnsets = [];
    JuOnsets = []; %Juice onsets
    SacSize = [];
    SacDur = [];
    %****** then allocate space and build structure
    SacX = zeros(Nstim,5);  % mat_time, ephy_time, stimulus
    StimK = 0;
    for zk = 1:size(Exp.D,1)
      if (length(Exp.D{zk}.PR.name) >= 6)  
        if (strcmp(Exp.D{zk}.PR.name(1:6),'Forage'))  % keep in mind, other trials could be image backgrounds
          if (Exp.D{zk}.PR.noisetype == 1) % Hartley
             %******************
             matstart = Exp.D{zk}.eyeData(1,6);  % register frame start clock
             matend = Exp.D{zk}.eyeData(end,6);  % register frame end clock
             eyeSmo = Exp.D{zk}.eyeSmo;    % vpx, smoothed data
             slist = Exp.D{zk}.slist;      % saccade times relative to data
             %***** loop through the saccade list for trial
             for k = 1:size(slist,1)
                 %***** get start and end
                 stasac = slist(k,4); % integer start time of saccade
                 endsac = slist(k,5); % integer end time of saccade
                 stx = eyeSmo(stasac,2);
                 sty = eyeSmo(stasac,3);
                 etx = eyeSmo(endsac,2); 
                 ety = eyeSmo(endsac,3);
                 %*********
                 ontime = matstart + slist(k,1); % time in seconds from trial start
                 offtime = matstart + slist(k,2); % time in seconds from trial start
                 %**********
                 %**** know when saccade starts relative to trial start
                 %**** so then determine what frame that is, and where
                 %**** the target is located, and if she went to it, and
                 %**** if it were pref or non-pref
                 %***********
                 z = find( (Exp.D{zk}.PR.ProbeHistory(:,4) >= ontime) & ...
                           (Exp.D{zk}.PR.ProbeHistory(:,4) < (ontime+MinSacDur)) );
                 if isempty(z)
                     disp('Error finding saccade in PR History, continue');
                     continue;
                 end
                 %**********
                 sacstart = z(1);  % frame when saccade is onset
                 SacX(StimK+sacstart,1) = -1;  %start by knowing its a saccade
                 SacX(StimK+sacstart,2) = (etx-stx); % x step
                 SacX(StimK+sacstart,3) = (ety-sty); % y step
                 %*************
                 z2 = find( (Exp.D{zk}.PR.ProbeHistory(:,4) >= offtime) & ...
                            (Exp.D{zk}.PR.ProbeHistory(:,4) < (offtime+MinSacDur)) );
                 if isempty(z2)
                     disp('Error finding saccade in PR History, continue');
                     continue;
                 end
                 %**********
                 sacend = z2(1);  % frame when saccade is onset
                 SacX(StimK+sacend,4) = -1;  %start by knowing its a saccad
                 %******* here you would determine if it goes to a target
                 if (Exp.D{zk}.PR.ProbeHistory(sacstart,3) > 0)  % a target was on
                     %******* do you go to the target
                     tx = Exp.D{zk}.PR.ProbeHistory(sacstart,1); % targ x pos
                     ty = Exp.D{zk}.PR.ProbeHistory(sacstart,2); % targ y pos
                     or = Exp.D{zk}.PR.ProbeHistory(sacstart,3); % targ ori, 1pref, 2 nonpref
                     %*******************
                     fdist1 = norm([stx-tx,sty-ty]);  % dist before saccade
                     fdist2 = norm([etx-tx,ety-ty]);  % dist after saccade
                     if (fdist1 > TargWin) && (fdist2 <= TargWin)
                         SacX(StimK+sacstart,1) = floor( or );
                         %******* determine if stimulus was off before
                         %******* saccade offset, only include those trials
                         OffBeforeEnd = 1;
                         if (size(Exp.D{zk}.PR.ProbeHistory,2) < 5)
                             JuOnsets = [JuOnsets ; NaN];
                         else
                             for ibk = sacend:-1:(sacend-3)
                                 if ( Exp.D{zk}.PR.ProbeHistory(sacend,5) == 1)
                                     OffBeforeEnd = 0;
                                 end
                             end
                             %****************
                             if (OffBeforeEnd)
                                 SacX(StimK+sacend,4) = floor( or );
                                 %******* next, if it was off, find re-onset
                                 ibk = sacend;
                                 while (ibk <= size(Exp.D{zk}.PR.ProbeHistory,1))
                                     if (Exp.D{zk}.PR.ProbeHistory(ibk,5) == 1)
                                         break;
                                     else
                                         ibk = ibk + 1;
                                     end
                                 end
                                 if (ibk < size(Exp.D{zk}.PR.ProbeHistory,1))
                                      or = Exp.D{zk}.PR.ProbeHistory(ibk,3); % targ ori, 1pref, 2 nonpref
                                      jtime = Exp.D{zk}.PR.ProbeHistory(ibk,4);
                                      if ((jtime-offtime) > 0.1)  % beyond normal jitter
                                          JuOnsets = [JuOnsets ; NaN];
                                      else
                                          SacX(StimK+ibk,5) = floor( or);  % post-sac stim ori
                                          JuOnsets = [JuOnsets ; jtime];
                                      end
                                 else
                                      JuOnsets = [JuOnsets ; NaN];
                                 end
                             else
                                 JuOnsets = [JuOnsets ; NaN]; 
                             end 
                         end
                         
                         %****** located stim re-onset on sac landing
%                          for ibk = 1:size(Exp.D{zk}.eyeData,1)
%                              if (Exp.D{zk}.eyeData(ibk,6) >= ontime)
%                                  break;
%                              end
%                          end
%                          while (ibk <= size(Exp.D{zk}.eyeData,1))
%                              if (Exp.D{zk}.eyeData(ibk,5) > 1) % state is 2 or 3
%                                  break;
%                              else
%                                  ibk = ibk + 1;
%                              end
%                          end
%                          if (ibk < size(Exp.D{zk}.eyeData,1))
%                             jtime = Exp.D{zk}.eyeData(ibk,6);
%                             %******* 
%                             if ( ((jtime-offtime) > 0.05) && ((jtime-offtime)<0.15) )
%                                 bzk = sacstart;
%                                 while (bzk <= size(Exp.D{zk}.PR.ProbeHistory,1))
%                                     if (Exp.D{zk}.PR.ProbeHistory(bzk,4) >= (jtime-0.09))
%                                         break;
%                                     else
%                                         bzk = bzk + 1;
%                                     end
%                                 end
%                                 if (bzk < size(Exp.D{zk}.PR.ProbeHistory,1))
%                                    SacX(StimK+bzk,5) = floor(or);
%                                 end
%                             end
%                             %************
%                             JuOnsets = [JuOnsets ; jtime];
%                          else
%                             JuOnsets = [JuOnsets ; NaN]; 
%                          end
%                          
                         
                         %******** identified a target saccade, when is
                         %******** sac onset compared to stim onset
                         SacOnsets = [SacOnsets ; ontime];
                         %****** working from that time, find target onset
                         izk = sacstart;
                         while (izk > 0)
                            if ~isnan(Exp.D{zk}.PR.ProbeHistory(izk,3)) 
                              izk = izk - 1;
                            else
                              break;
                            end
                         end
                         if (izk > 0)
                             stimtime = Exp.D{zk}.PR.ProbeHistory(izk,4);
                             StimOnsets = [StimOnsets ; stimtime];
                             SacSize = [SacSize ; norm([(etx-stx),(ety-sty)])];
                             SacOffsets = [SacOffsets ; offtime];
                             SacDur = [SacDur ; (offtime-ontime)];
                         else
                             StimOnsets = [StimOnsets ; NaN];
                             SacSize = [SacSize ; NaN];
                             SacOffsets = [SacOffsets ; NaN];
                             SacDur = [SacDur ; NaN];
                         end     
                         %***************************************
                     end
                 end
             end
             %******** build stimulus matrix, step to next trial
             StimK = StimK + size(Exp.D{zk}.PR.NoiseHistory,1);
             %*************************
          end
        end
      end
      if (mod(zk,10)==0)
         disp(sprintf('Processing trial %d for noise history',zk));
      end
    end
    SacX = SacX(1:StimK,:);
    %*********************************
    
    
    %******** plot the reactions times in the task
    H = figure;
    set(H,'position',[100 100 800 800]);
    subplot('position',[0.1 0.4 0.4 0.5]);
    plot((SacOnsets - StimOnsets),1:size(StimOnsets,1),'k.'); hold on;
    xlabel('Reaction time');
    ylabel('Trial');
    axis([0 1 0 size(StimOnsets,1)]);
    title('Plot of target reaction times in Forage 3');
    nanmean( (SacOnsets - StimOnsets))
    nanstd( (SacOnsets - StimOnsets))
    nanmedian( (SacOnsets - StimOnsets) )
    z = find( (SacOnsets - StimOnsets) < 0.1 );  % frac sacs under 100ms
    length(z)/length(SacOnsets)
    %********
    subplot('position',[0.6 0.7 0.3 0.25]);
    ivx = 0:0.025:0.8;
    yvx = hist((SacOnsets-StimOnsets),ivx);
    bar(ivx,yvx);
    xlabel('Reaction time');
    ylabel('Sac counts');
    %*****************
    subplot('position',[0.6 0.4 0.3 0.25]);
    z = find( (SacSize > 1) & ((SacOnsets-StimOnsets) < 0.4) );
    plot(SacSize(z),(SacOnsets(z)-StimOnsets(z)),'k.');
    xlabel('Sac Size');
    ylabel('Sac RT');
    %*****************
    subplot('position',[0.6 0.1 0.3 0.25]);
    ivx = 0:0.005:0.10;
    yvx = hist(SacDur,ivx);
    bar(ivx,yvx);
    xlabel('Sac Dur');
    ylabel('Sac counts');
    %***************
    subplot('position',[0.1 0.05 0.3 0.3]);
    z = find(SacDur < 0.1);
    plot(SacDur(z),SacSize(z),'k.');
    xlabel('Sac Dur');
    ylabel('Sac Size');
    %*******************************************
    
    figure;
    plot((JuOnsets-SacOffsets),1:size(StimOnsets,1),'k.');
    xlabel('Reward latency from sac offset');
    ylabel('Trial');
    
return;
