function [RFsets,ztXX,Zx,Zy,mino,maxo] = PlotForageSpatialKernel_RF(StimX,StimY,GRID)
% function PlotForageSpatialKernel(Exp)
%   input: takes StimX and StimY (stim history and spike counts)
%     StimX: [Nframes,3 + NT] - fields, matlab time, ephys time, stim, (if stim)
%     StimY: [Nframes,1] - spike counts per video frame
%
%******** use reverse correlation of stim vs no stim to find temporal
%******** the temporal kernel from the data stream
SampRate = median(diff(StimX(:,1)));  % median screen flip time
DTA = -5;   % stim frames, relative to stim onset, must be negative
DTB = 18;   % stim frames after onset, must be positive 
XX = DTA:DTB;
tXX = XX * (SampRate * 1000);  % put into ms
SN = size(StimX,1);
%**** transform some parameters of the grid
TT = size(XX,2);
Nx = floor(GRID.box(3)/GRID.div);
Ny = floor(GRID.box(4)/GRID.div);
NT = (Nx * Ny);
BaseX = GRID.box(1) - (GRID.box(3)/2);
BaseY = GRID.box(2) - (GRID.box(4)/2);
Div = GRID.div;
Zx = BaseX + (1:Nx)*Div;
Zy = BaseY + (1:Ny)*Div;
%***************
mcounts = zeros(NT,TT);   %response to stimulus
mcounts2 = zeros(NT,TT);
counts = zeros(NT,1);
msem = zeros(NT,TT);
for k = 1:SN
    svec = StimX(k,3:end);
    z = find( svec ~= 0);   
    if ~isempty(z)
        if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnts = StimY((k+XX)',1)';
             for it = z
               mcounts(it,:) = mcounts(it,:) + spcnts;
               mcounts2(it,:) = mcounts2(it,:) + spcnts .^ 2;
               counts(it) = counts(it) + 1;
             end
        end
    end
    if (mod(k,floor(SN/10)) == 0)
       disp(sprintf('Processing %5.1f percent comp',(k*100/SN))); 
    end
end
for it = 1:NT
    if counts(it) > 0
       mcounts(it,:) = mcounts(it,:) / counts(it);  % mean
       msem(it,:) = sqrt(  ((mcounts2(it,:)/counts(it)) - mcounts(it,:).^2) / counts(it));  % 1 sem 
       mcounts(it,:) = mcounts(it,:) / SampRate;
       msem(it,:) = msem(it,:) / SampRate;
    else
       mcounts = NaN(size(mcounts));
       msem = NaN(size(mcounts));
    end
end
%******* if no counts there is a problem, return
if isnan(mcounts) 
   disp('Failed to estimate any temporal kernel, stopping'); 
   return;  % otherwise return, if no temporal kernel there is a problem
end

%********* for starts, store and return to plot in other routine
KN = 7;
RFsets = cell(1,KN);
mino = min(min(mcounts));
maxo = max(max(mcounts));
for it = 1:(KN-1)
   ito = (5-DTA)+it;
   svec = flipud( reshape(squeeze(mcounts(:,ito)),Nx,Ny)' );
   RFsets{it} = svec;
   ztXX(it) = -tXX(ito);
end
%********************
   
return;
