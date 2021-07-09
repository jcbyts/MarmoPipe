function PlotForageSpatialKernel(StimX,StimY,GRID,FileTag)
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
    %z = find( svec > 0);
    %z = find( svec < 0);
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
%********* for starts, just plot entire correlation for all frames
hf = figure;
set(hf,'position',[100 100 800 800]);
KN = ceil( sqrt(TT) );
mino = min(min(mcounts));
maxo = max(max(mcounts));
for it = 1:TT
   subplot(KN,KN,it);
   svec = flipud( reshape(squeeze(mcounts(:,it)),Nx,Ny)' );
   imagesc(Zx,Zy,svec,[mino maxo]); hold on;
   plot([Zx(1),Zx(end)],[0,0],'k-');
   plot([0,0],[Zy(1),Zy(end)],'k-');  
   title(sprintf('%4.1f ms',tXX(it)));
end
%******* plot the counts per bin
maxo = max(counts);
mino = min(counts);
[mino maxo]
subplot(KN,KN,(KN*KN));
svec = flipud( reshape(counts,Nx,Ny)');
imagesc(Zx,Zy,svec,[0 maxo]); hold on;
plot([Zx(1),Zx(end)],[0,0],'k-');
plot([0,0],[Zy(1),Zy(end)],'k-');
title('Counts');

%********* for starts, just plot entire correlation for all frames
hf2 = figure;
set(hf2,'position',[300 300 800 100]);
KN = 9; %ceil( sqrt(TT) );
mino = min(min(mcounts));
maxo = max(max(mcounts));
for it = 1:(KN-1)
   subplot(1,KN,(KN-it));
   ito = (3-DTA)+it;
   svec = flipud( reshape(squeeze(mcounts(:,ito)),Nx,Ny)' );
   imagesc(Zx,Zy,svec,[mino maxo]); hold on;
   plot([-11,11],[0,0],'k-');
   plot([0,0],[-11,11],'k-');
   axis off;
   h = title(sprintf('%4.1f',-tXX(ito)));
end
subplot(1,KN,KN);
imagesc(Zx,Zy,ones(size(svec))*mino,[mino maxo]); hold on;
axis off;
h = colorbar;
%******* plot the counts per bin

   
return;
