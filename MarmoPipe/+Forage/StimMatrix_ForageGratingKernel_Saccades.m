function [SacX] = StimMatrix_ForageGratingKernel_Saccades(Exp,SPclust)
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
%*************

    sp = Exp.sp{SPclust};
    
    %**** go through trials, building a stimulus representations
    %**** in the coordinate frame of spike times ***************
    %**** at this time, no consideration of saccades or probe
    %**** stimulus timing is considered (later)
    %***
    %***** first scroll trials to compute size of stim array in time
    Nstim = 0;
    for zk = 1:size(Exp.D,1)
       if strcmp(Exp.D{zk}.PR.name,'Forage')
          if (Exp.D{zk}.PR.noisetype == 1) % Hartley
            Nstim = Nstim + size(Exp.D{zk}.PR.NoiseHistory,1);
          end
       end
    end
     
    %****** then allocate space and build structure
    SacX = zeros(Nstim,3);  % mat_time, ephy_time, stimulus
    StimK = 0;
    for zk = 1:size(Exp.D,1)
       if (strcmp(Exp.D{zk}.PR.name,'Forage'))  % keep in mind, other trials could be image backgrounds
          if (Exp.D{zk}.PR.noisetype == 1) % Hartley
             %******************
              
             %******************
             matstart = Exp.D{zk}.eyeData(1,6);  % register frame start clock
             matend = Exp.D{zk}.eyeData(end,6);  % register frame end clock
             epystart = Exp.D{zk}.START_EPHYS;   % matching strobe time ephys start
             epyend = Exp.D{zk}.END_EPHYS;       % matching strobe time ephys end
             %******** build stimulus matrix in spike time coordinates
             timvec = Exp.D{zk}.PR.NoiseHistory(:,1);  % delayed by 8 ms, fix in stim
             stimvec = Exp.D{zk}.PR.NoiseHistory(:,2);
             epyvec = (timvec-matstart)/(matend-matstart);
             epyvec = epystart + epyvec * (epyend - epystart);
             %***************************
             ia = StimK+1;
             ib = StimK+size(timvec,1);
             StimX(ia:ib,1) = timvec;
             StimX(ia:ib,2) = epyvec;
             StimX(ia:ib,3) = stimvec;
             StimX(ia:ib,4) = (stimvec > 0);
             %****** now compute counts via frame
             for it = 2:size(epyvec,1)
                  zz = find( (sp.st >= epyvec(it-1)) & (sp.st < epyvec(it) ) );
                  StimY(ia+(it-1),1) = size(zz,1);
             end
             StimK = StimK + size(timvec,1);
             %*************************
          end
       end
       if (mod(zk,10)==0)
         disp(sprintf('Processing trial %d for noise history',zk));
       end
    end
    StimX = StimX(1:StimK,:);
    StimY = StimY(1:StimK,:);
    %*********************************
    
return;
