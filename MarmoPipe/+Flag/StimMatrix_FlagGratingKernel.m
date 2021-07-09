function [StimX,StimY] = StimMatrix_ForageGratingKernel(Exp)
% function StimMatrix_ForageGratingKernel(Exp)
%   input: takes Exp struct as the input
%  Exp fields:
%    D: {100x1 cell}  - data per trial
%    S: [1x1 struct]  - rig settings struct
%   sp: [1x1 struct]  - spikes from spike sorting
%   cp: (will be a 60hz filtered LFP at 1000 hz)
%   ey: (will be a struct for eye position)
% Outputs:
%   StimX: [Nframes,4] - fields, matlab time, ephys time, stim, (if stim)
%   StimY: [Nframes,1] - spike counts per video frame
%*************

    %**** go through trials, building a stimulus representations
    %**** in the coordinate frame of spike times ***************
    %**** at this time, no consideration of saccades or probe
    %**** stimulus timing is considered (later)
    %***
    %***** first scroll trials to compute size of stim array in time
    Nstim = 0;
    for zk = 1:size(Exp.D,1)
       if strcmp(Exp.D{zk}.PR.name,'Flag')
          Nstim = Nstim + size(Exp.D{zk}.PR.NoiseHistory,1);
       end
    end
    %****** then allocate space and build structure
    StimX = zeros(Nstim,5);  % mat_time, ephy_time, stimulus
    StimY = zeros(Nstim,1);  % spike counts per frame epoch
    StimK = 0;
    for zk = 1:size(Exp.D,1)
       if strcmp(Exp.D{zk}.PR.name,'Flag') 
         %******************
         matstart = Exp.D{zk}.eyeData(1,1);  % register frame start clock
         matend = Exp.D{zk}.eyeData(end,6);  % register frame end clock
         epystart = Exp.D{zk}.START_EPHYS;   % matching strobe time ephys start
         epyend = Exp.D{zk}.END_EPHYS;       % matching strobe time ephys end
         %******** identify target orientation on trial and when target is
         %******** visible to possibly influence the firing ***********
         onstim = NaN;
         itargori = NaN;
         if (Exp.D{zk}.PR.error == 0)  % only correct trials
            z = find( Exp.D{zk}.eyeData(:,5) == 2);
            if ~isempty(z)
              onstim = Exp.D{zk}.eyeData(z(1),1);
            end
            %**** 
            z = find( Exp.D{zk}.eyeData(:,5) >= 5);
            if ~isempty(z)
              offstim = Exp.D{zk}.eyeData(z(1),1);
            end
            %****
            OriNum = Exp.D{zk}.P.orinum; % number of orientations for targets
            OriDiff = (180/OriNum);  % angle diff between orientations
            itargori = 1+floor(Exp.D{zk}.PR.targori/OriDiff);  % assumes 12 orientations
            %******
         end
         if ~isnan(onstim)
            onstim = (onstim-matstart)/(matend-matstart);
            onstim = epystart + onstim * (epyend - epystart);
            offstim = (offstim-matstart)/(matend-matstart);
            offstim = epystart + offstim * (epyend - epystart);
         end        
         %******** build stimulus matrix in spike time coordinates
         if ~isempty(Exp.D{zk}.PR.NoiseHistory)  % could be empty if she aborts trial
             timvec = Exp.D{zk}.PR.NoiseHistory(:,1);  % delayed by 8 ms, fix in stim
             stimvec = Exp.D{zk}.PR.NoiseHistory(:,2);
             epyvec = (timvec-matstart)/(matend-matstart);
             epyvec = epystart + epyvec * (epyend - epystart);
             %***************************
             ia = StimK+1;
             ib = StimK+size(timvec,1);
             snip = ia:ib;
             StimX(snip,1) = timvec;
             StimX(snip,2) = epyvec;
             StimX(snip,3) = stimvec;
             StimX(snip,4) = (stimvec > 0);
             %****** mark any frames where target stim was present, record
             %****** its orientation on those frames
             if ~isnan(onstim) && ~isnan(itargori)
                z = find( (epyvec > onstim) & (epyvec < offstim) );
                if ~isempty(z)
                    StimX(snip(z),5) = itargori;
                end
             end
             %****** now compute counts via frame
             for it = 2:size(epyvec,1)
                  zz = find( (Exp.sp.st >= epyvec(it-1)) & (Exp.sp.st < epyvec(it) ) );
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
