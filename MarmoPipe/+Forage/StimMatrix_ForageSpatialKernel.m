function [StimX,StimXP,StimXN,StimY] = StimMatrix_ForageSpatialKernel(Exp,GRID,SPclust)
% function StimMatrix_ForageSpatialKernel(Exp)
%   input: takes Exp struct as the input
%          GRID with fields:
%             GRID.box = [cent x, cent y, width, height]
%             GRID.div = size of spatial pixel
%          SPClust - number of spike cluster to analyze
%
%  Exp fields:
%    D: {100x1 cell}  - data per trial
%    S: [1x1 struct]  - rig settings struct
%   sp: [1x1 struct]  - spikes from spike sorting
%   cp: (will be a 60hz filtered LFP at 1000 hz)
%   ey: (will be a struct for eye position)
% Outputs:
%            note, StimX,StimY are relative to the eye position (eye corrected)
%   StimX: [Nframes,3+NK] - fields, matlab time, ephys time, NK = size of grid
%   StimXP:  identical to StimX, but only white flashes
%   StimXN:  identical to StimY, but only black flashes
%   StimY: [Nframes,1] - spike counts per video frame
%*************

    sp = Exp.sp{SPclust};
    StimX = [];
    StimXP = [];
    StimXN = [];
    StimY = [];
    
    %**** transform some parameters of the grid
    Nx = floor(GRID.box(3)/GRID.div);
    Ny = floor(GRID.box(4)/GRID.div);
    NT = (Nx * Ny);
    BaseX = GRID.box(1) - (GRID.box(3)/2);
    BaseY = GRID.box(2) - (GRID.box(4)/2);
    Div = GRID.div;
    %***********************
    
    %**** go through trials, building a stimulus representations
    %**** in the coordinate frame of spike times ***************
    %**** at this time, no consideration of saccades or probe
    %**** stimulus timing is considered (later)
    %***
    %***** first scroll trials to compute size of stim array in time
    Nstim = 0;
    disp('Testing for trials with spatial mapping');
    for zk = 1:size(Exp.D,1)
      if (length(Exp.D{zk}.PR.name) >= 6)
        if strcmp(Exp.D{zk}.PR.name(1:6),'Forage')
          if (Exp.D{zk}.PR.noisetype == 2) % Spatial
            disp(sprintf('Including trial %d',zk));    
            Nstim = Nstim + size(Exp.D{zk}.PR.NoiseHistory,1);
          end
        end
      end
    end
    disp(sprintf('Collected %d frames',Nstim));
    if (Nstim == 0)
        return;
    end
    
    %****** then allocate space and build structure
    StimX = zeros(Nstim,3+NT-1);  % mat_time, ephy_time, stimulus
    StimXP = zeros(Nstim,3+NT-1);  % mat_time, ephy_time, stimulus
    StimXN = zeros(Nstim,3+NT-1);  % mat_time, ephy_time, stimulus
    StimY = zeros(Nstim,1);  % spike counts per frame epoch
    StimK = 0;
    for zk = 1:size(Exp.D,1)
      if (length(Exp.D{zk}.PR.name) >= 6)
        if (strcmp(Exp.D{zk}.PR.name(1:6),'Forage'))  % keep in mind, other trials could be image backgrounds
          if (Exp.D{zk}.PR.noisetype == 2) % Spatial
             %******************
             matstart = Exp.D{zk}.eyeData(1,6);  % register frame start clock
             matend = Exp.D{zk}.eyeData(end,6);  % register frame end clock
             epystart = Exp.D{zk}.START_EPHYS;   % matching strobe time ephys start
             epyend = Exp.D{zk}.END_EPHYS;       % matching strobe time ephys end
             eyesmo = Exp.D{zk}.eyeSmo;          % recorded eye position
             %******** build stimulus matrix in spike time coordinates
             timvec = Exp.D{zk}.PR.NoiseHistory(:,1);  % delayed by 8 ms, fix in stim
             epyvec = (timvec-matstart)/(matend-matstart);
             epyvec = epystart + epyvec * (epyend - epystart);
             %***************************
             ia = StimK+1;
             ib = StimK+size(timvec,1);
             StimX(ia:ib,1) = timvec;
             StimX(ia:ib,2) = epyvec;
             StimXP(ia:ib,1:2) = StimX(ia:ib,1:2);
             StimXN(ia:ib,1:2) = StimX(ia:ib,1:2);
             %******* loop through to compute spatial positions
             KT = size(timvec,1);
             svec = zeros(KT,NT);
             svecp = zeros(KT,NT);
             svecn = zeros(KT,NT);
             %****
             for k = 1:KT
                 framet = timvec(k) - matstart;
                 z = find( eyesmo(:,1) >= framet );
                 if ~isempty(z)
                    ex = eyesmo(z(1),2);
                    ey = eyesmo(z(1),3);
                    for sk = 1:Exp.D{zk}.PR.noiseNum
                       sx = Exp.D{zk}.PR.NoiseHistory(k,2+((sk-1)*2));  % x pos
                       sy = Exp.D{zk}.PR.NoiseHistory(k,3+((sk-1)*2));  % y pos
                       cval = ((mod(sk,2)*2) - 1);  % white is 1, black -1
                       %********* determine grid pixel and mark with val
                       sxx = sx - ex;  % recentered by eye pos
                       syy = sy - ey;
                       %****
                       ix = 1 + floor((sxx - BaseX -(Div/2))/Div);
                       iy = 1 + floor((syy - BaseY -(Div/2))/Div);
                       if ( (ix >= 1) && (ix <= Nx) && (iy >= 1) && (iy <= Ny) )
                          vox = (iy-1)*Nx + ix;
                          svec(k,vox) = svec(k,vox) + cval;
                          if (cval > 0)
                             svecp(k,vox) = svecp(k,vox) + cval;
                          else
                             svecn(k,vox) = svecn(k,vox) + cval;
                          end
                       end
                       %*******************************************
                    end
                 end
             end
             %****
             StimX(ia:ib,3:(NT+2)) = svec;
             StimXP(ia:ib,3:(NT+2)) = svecp;
             StimXN(ia:ib,3:(NT+2)) = svecn;
             %******************************************
             %****** now compute counts via frame
             for it = 2:size(epyvec,1)
                  zz = find( (sp.st >= epyvec(it-1)) & (sp.st < epyvec(it) ) );
                  StimY(ia+(it-1),1) = size(zz,1);
             end
             StimK = StimK + size(timvec,1);
             %*************************
          end
        end
      end
      if (mod(zk,10)==0)
         disp(sprintf('Processing trial %d for noise history',zk));
      end
    end
    StimX = StimX(1:StimK,:);
    StimXP = StimXP(1:StimK,:);
    StimXN = StimXN(1:StimK,:);
    StimY = StimY(1:StimK,:);
    disp(sprintf('Processing finished'));
    %*********************************

return;
