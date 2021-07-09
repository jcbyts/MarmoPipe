function PlotForageSaccadeKernel(SacX,StimX,StimY,FileTag,SO)
% function PlotForageSaccadeKernel(SacX,StimX,StimY,FileTag)
%   input: takes SacX and StimY (sac history and spike counts)
%     SacX:  [Nframes,3] - fields, for saccade timing
%     StimX: [Nframes,4] - fields, matlab time, ephys time, stim, (if stim)
%     StimY: [Nframes,1] - spike counts per video frame
%     FileTag: name of file
%     SO:  saccade onset lock - 1 - sac onset, 4 -sac offset, 
%                               5 - stim reappear inside RF after sac
%   output: plot of the temporal and tuning kernels

%******** use reverse correlation of stim vs no stim to find temporal
%******** the temporal kernel from the data stream
SampRate = median(diff(StimX(:,1)));  % median screen flip time
DTA = -12;   % stim frames, relative to stim onset, must be negative
DTB = 24;   % stim frames after onset, must be positive 
XX = DTA:DTB;
tXX = XX * (SampRate * 1000);  % put into ms
SN = size(StimX,1);
%*******
STARNG = -17; % frames prior to sac onset
ENDRNG = -2; %10;
REVTEMPK = false;
PrefSacs = [11:12,1]; %7:12]; 
NonPrefSacs = [5:7]; 
SameTarg = 2;
%***************
mcounts = zeros(size(XX));   %response to stimulus
mcounts2 = zeros(size(XX));
counts = 0;
for k = 1:SN
    if (SacX(k,SO))  %only analyze inside range, use saccade onsets instead
        if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnts = StimY((k+XX)',1)';
             mcounts = mcounts + spcnts;
             mcounts2 = mcounts2 + spcnts .^ 2;
             counts = counts + 1;
        end
    end
end
if counts > 0
   mcounts = mcounts / counts;  % mean
   msem = sqrt(  ((mcounts2/counts) - mcounts.^2) / counts);  % 1 sem 
   mcounts = mcounts / SampRate;
   msem = msem / SampRate;
else
   mcounts = NaN(size(mcounts));
   msem = NaN(size(mcounts));
end
%***************
Amcounts = zeros(size(XX));   %response to stimulus
Amcounts2 = zeros(size(XX));
Acounts = 0;
for k = 1:SN
    if (ismember(SacX(k,SO),PrefSacs))  %only analyze inside range, use saccade onsets instead
      %*****
      dipo = 1;
      if (SameTarg > 0)
         if (SameTarg == 1)
             if (SacX(k,1) == SacX(k,5)) % same target
                 dipo = 1;
             else
                 dipo = 0;
             end
         else
             if (SacX(k,1) == SacX(k,5)) % same target
                 dipo = 0;
             else
                 dipo = 1;
             end    
         end
      end
      %*******
      if (dipo)
        if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnts = StimY((k+XX)',1)';
             Amcounts = Amcounts + spcnts;
             Amcounts2 = Amcounts2 + spcnts .^ 2;
             Acounts = Acounts + 1;
        end
      end
      %*********
    end
end
if Acounts > 0
   Amcounts = Amcounts / Acounts;  % mean
   Amsem = sqrt(  ((Amcounts2/Acounts) - Amcounts.^2) / Acounts);  % 1 sem 
   Amcounts = Amcounts / SampRate;
   Amsem = Amsem / SampRate;
else
   Amcounts = NaN(size(Amcounts));
   Amsem = NaN(size(Amcounts));
end
%***************
Bmcounts = zeros(size(XX));   %response to stimulus
Bmcounts2 = zeros(size(XX));
Bcounts = 0;
for k = 1:SN
    if (ismember(SacX(k,SO),NonPrefSacs))  %only analyze inside range, use saccade onsets instead
      %*****
      dipo = 1;
      if (SameTarg > 0)
         if (SameTarg == 1)
             if (SacX(k,1) == SacX(k,5)) % same target
                 dipo = 1;
             else
                 dipo = 0;
             end
         else
             if (SacX(k,1) == SacX(k,5)) % same target
                 dipo = 0;
             else
                 dipo = 1;
             end    
         end
      end
      %*******
      if (dipo)
        if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnts = StimY((k+XX)',1)';
             Bmcounts = Bmcounts + spcnts;
             Bmcounts2 = Bmcounts2 + spcnts .^ 2;
             Bcounts = Bcounts + 1;
        end
      end
      %*********
    end
end
if Bcounts > 0
   Bmcounts = Bmcounts / Bcounts;  % mean
   Bmsem = sqrt(  ((Bmcounts2/Bcounts) - Bmcounts.^2) / Bcounts);  % 1 sem 
   Bmcounts = Bmcounts / SampRate;
   Bmsem = Bmsem / SampRate;
else
   Bmcounts = NaN(size(Bmcounts));
   Bmsem = NaN(size(Bmcounts));
end
%********** compute the null distribution from blanks
pcounts = zeros(size(XX));   %response to blanks
pcounts2 = zeros(size(XX));
bcounts = 0;
for k = 1:SN
    if (SacX(k,SO)==0)  % get response to blanks
        if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnts = StimY((k+XX)',1)';
             pcounts = pcounts + spcnts;
             pcounts2 = pcounts2 + spcnts .^ 2;
             bcounts = bcounts + 1;
        end
    end
end
if bcounts > 0
   pcounts = pcounts / bcounts;  % mean
   psem = sqrt(  ((pcounts2/bcounts) - pcounts.^2) / bcounts);  % 1 sem 
   pcounts = pcounts / SampRate;
   psem = psem / SampRate;
else
   pcounts = NaN(size(pcounts));
   psem = NaN(size(pcounts));
end


%******* if no counts there is a problem, return
if isnan(mcounts(1)) || isnan(pcounts(1))
   disp('Failed to estimate any saccade temporal kernel, stopping'); 
   return;  % otherwise return, if no temporal kernel there is a problem
else
   disp('Finished saccade trig avg');
end

%******* compute summary stats ******
%********* then plot results ********
hf = figure;
set(hf,'position',[100 80 900 900]);
%**** plot the sac triggered temporal kernel for now
subplot('position',[0.1 0.65 0.4 0.25]);
ymean = mean([mcounts(1:abs(DTA)),pcounts(1:abs(DTA))]);
hh = plot(tXX,Amcounts,'r-'); hold on;
set(hh,'LineWidth',2);
plot(tXX,(Amcounts+(2*Amsem)),'r:');
plot(tXX,(Amcounts-(2*Amsem)),'r:');
hh = plot(tXX,Bmcounts,'b-'); hold on;
set(hh,'LineWidth',2);
plot(tXX,(Bmcounts+(2*Bmsem)),'b:');
plot(tXX,(Bmcounts-(2*Bmsem)),'b:');
plot([tXX(1),tXX(end)],[ymean,ymean],'k-');
%******* replot the other saccades (on avg smaller ones)
hh = plot(tXX,mcounts,'k-'); hold on;
set(hh,'LineWidth',2);
plot(tXX,(mcounts+(2*msem)),'k:');
plot(tXX,(mcounts-(2*msem)),'k:');
plot([tXX(1),tXX(end)],[ymean,ymean],'k-');
%****************
axis tight;
V = axis;
plot([0,0],[V(3),V(4)],'k-');
xlabel('Time (ms)');
ylabel('Rate (sp/s)');
title(sprintf('Sac-Trig %s: (%d)(%d,%d)',FileTag,counts,Acounts,Bcounts));
   
%******** Run stim-triggered analysis to get temporal kernel
%******** use reverse correlation of stim vs no stim to find temporal
%******** the temporal kernel from the data stream
SampRate = median(diff(StimX(:,1)));  % median screen flip time
DTA = -5;   % stim frames, relative to stim onset, must be negative
DTB = 20;   % stim frames after onset, must be positive 
XX = DTA:DTB;
tXX = XX * (SampRate * 1000);  % put into ms
SN = size(StimX,1);
%***************
mcounts = zeros(size(XX));   %response to stimulus
mcounts2 = zeros(size(XX));
counts = 0;
for k = 1:SN
    if (StimX(k,4))  %only analyze inside range
        if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnts = StimY((k+XX)',1)';
             mcounts = mcounts + spcnts;
             mcounts2 = mcounts2 + spcnts .^ 2;
             counts = counts + 1;
        end
    end
end
if counts > 0
   mcounts = mcounts / counts;  % mean
   msem = sqrt(  ((mcounts2/counts) - mcounts.^2) / counts);  % 1 sem 
   mcounts = mcounts / SampRate;
   msem = msem / SampRate;
else
   mcounts = NaN(size(mcounts));
   msem = NaN(size(mcounts));
end
%********** compute the null distribution from blanks
pcounts = zeros(size(XX));   %response to blanks
pcounts2 = zeros(size(XX));
bcounts = 0;
for k = 1:SN
    if (~StimX(k,4))  % get response to blanks
        if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnts = StimY((k+XX)',1)';
             pcounts = pcounts + spcnts;
             pcounts2 = pcounts2 + spcnts .^ 2;
             bcounts = bcounts + 1;
        end
    end
end
if bcounts > 0
   pcounts = pcounts / bcounts;  % mean
   psem = sqrt(  ((pcounts2/bcounts) - pcounts.^2) / bcounts);  % 1 sem 
   pcounts = pcounts / SampRate;
   psem = psem / SampRate;
else
   pcounts = NaN(size(pcounts));
   psem = NaN(size(pcounts));
end
%******* if no counts there is a problem, return
if isnan(mcounts(1)) || isnan(pcounts(1))
   disp('Failed to estimate any visual temporal kernel, stopping'); 
   return;  % otherwise return, if no temporal kernel there is a problem
else
    disp('Finished visual temporal kernel');
end

%******* Find what you consider a significant temporal kernel
trepa = 0;
trepb = 0;
tcrit = 8;  %threshold for significance, in sem
if (~REVTEMPK)
    k = 1;
    while k <= size(XX,2)
        if (trepa == 0)  % look for first sig time point
            if ((mcounts(k)-pcounts(k)) > tcrit*psem(k))
                trepa = k;  % flag first sig time point
            end
        else
            if (trepb == 0)  % find first time point after trepa not sig
                if ((mcounts(k)-pcounts(k)) < tcrit*psem(k))
                    trepb = k;
                end
            end
        end
        k = k + 1;
    end
else
    k = size(XX,2);
    while k >= 1
        if (trepb == 0)  % look for first sig time point
            if ((mcounts(k)-pcounts(k)) > tcrit*psem(k))
                trepb = k;  % flag first sig time point
            end
        else
            if (trepa == 0)  % find first time point after trepa not sig
                if ((mcounts(k)-pcounts(k)) < tcrit*psem(k))
                    trepa = k;
                end
            end
        end
        k = k - 1;
    end
end
%*******
if (trepa > 0) && (trepb == 0)
    trepb = size(XX,2);  % stay sig long term
end

%******* compute summary stats ******
%********* then plot results ********
%**** plot the temporal kernel for now
subplot('position',[0.6 0.65 0.35 0.25]);
ymean = mean([mcounts(1:abs(DTA)),pcounts(1:abs(DTA))]);
plot(tXX,mcounts,'b-'); hold on;
plot(tXX,(mcounts+(2*msem)),'b--');
plot(tXX,(mcounts-(2*msem)),'b--');
plot(tXX,(pcounts+(tcrit*psem)),'b:');
plot(tXX,(pcounts-(tcrit*psem)),'b:');
plot([tXX(1),tXX(end)],[ymean,ymean],'k-');
axis tight;
V = axis;
plot([0,0],[V(3),V(4)],'k-');
plot([tXX(trepa),tXX(trepa)],[V(3),V(4)],'k--');
plot([tXX(trepb),tXX(trepb)],[V(3),V(4)],'k--');
xlabel('Time (ms)');
ylabel('Rate (sp/s)');
title(sprintf('Temporal Kernel: %s',FileTag));

%********* if no significant temporal kernel then return ***
if ~trepa || ~trepb
    disp('No significant temporal kernel found, stopping');
    return;
end

%****** Now using that sig temporal kernel, look for feature selectivity
%****** across spat freq and orientation *********************
StimList = unique(StimX(:,3));  %72 or 32, but plus 1
if (length(StimList) == 73)
  SpatOris = (0:11) * 15;   % orientations, figure way to get auto from Exp
  SpatFrqs = exp(log(0.5):((log(16)-log(0.5))/5):log(16));
  MyTick = [0.5:((16-0.5)/5):16];
  MyTickLabel = {'0.5','1','2','4','8','16'};
else
  if (length(StimList) == 33)   
      SpatOris = (0:7) * 22.5;   % orientations, figure way to get auto from Exp
      SpatFrqs = exp(log(2):((log(16)-log(2))/3):log(16));
      MyTick = [2:((16-2)/3):16];
      MyTickLabel = {'2','4','8','16'}; 
  else
     if (length(StimList) == 49)
       SpatOris = (0:11) * 15;   % orientations, figure way to get auto from Exp
       SpatFrqs = exp(log(2):((log(16)-log(2))/3):log(16));
       MyTick = [2:((16-2)/3):16];
       MyTickLabel = {'2','4','8','16'}; 
     else
      if (length(StimList) == 61)
           SpatOris = (0:11) * 15;   % orientations, figure way to get auto from Exp
           SpatFrqs = exp(log(1):((log(16)-log(1))/3):log(16));
           MyTick = [1:((16-1)/4):16];
           MyTickLabel = {'1','2','4','8','16'}; 
       else
            disp('size mismatch, return');
            return;
       end
     end
  end
end
%******
NOri = size(SpatOris,2);
NSpf = size(SpatFrqs,2);
NTot = NOri * NSpf;
%**** run same kind of averages, but broken out by stim type
kcounts = zeros(1,NTot);
stcounts = zeros(1,NTot);
stcounts2 = zeros(1,NTot);
for k = 1:SN
    if (StimX(k,4))  %only analyze inside range
        stimo = StimX(k,3); % stimulus type
        if ( (stimo>0) && (stimo<=NTot))
           if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnt = mean( StimY((k+XX(trepa:trepb))',1)');  % mean over epoch
             stcounts(stimo) = stcounts(stimo) + spcnt;
             stcounts2(stimo) = stcounts2(stimo) + (spcnt ^ 2);
             kcounts(stimo) = kcounts(stimo) + 1;
           end
        end
    end
end
if kcounts > 0
   stcounts = stcounts ./ kcounts;  % mean
   stsem = sqrt(  ((stcounts2 ./ kcounts) - stcounts.^2) ./ kcounts);  % 1 sem 
   stcounts = stcounts / SampRate;
   stsem = stsem / SampRate;
else
   stcounts = NaN(size(stcounts));
   stsem = NaN(size(stcounts));
end

%******* plot out the stimulus selectivity over sig temporal epoch
subplot('position',[0.1 0.35 0.25 0.20]);
uu = reshape(stcounts,[NOri NSpf]);
su = reshape(stsem,[NOri NSpf]);
maxo = max(max( uu + su ));
mino = min(min( uu - su ));
mino = mino * 0.95;
imagesc(SpatOris,SpatFrqs,uu',[mino maxo]); colormap('gray'); h = colorbar;
xlabel('Orientation');
ylabel('Spat Freq');
set(gca,'YTick',MyTick,'YTickLabel',MyTickLabel);
title('Rate(sp/s)');

%***** plot orientation average
subplot('position',[0.4 0.35 0.225 0.20]);
otune = mean(uu');
otune = [otune otune];
sotune = mean(su');
sotune = [sotune sotune];
xo = [SpatOris (SpatOris+180)];
hh = plot(xo,otune,'k-'); hold on;
set(hh,'Linewidth',2);
plot(xo,otune + (2*sotune),'k:'); hold on;
plot(xo,otune - (2*sotune),'k:'); hold on;
plot(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('Orientation (degs)');
ylabel('Rate (sp/s)');

%***** plot spf average
subplot('position',[0.7 0.35 0.225 0.20]);
stune = mean(uu);
sstune = mean(su);
xo = SpatFrqs; 
hh = semilogx(xo,stune,'k-'); hold on;
set(hh,'Linewidth',2);
semilogx(xo,stune + (2*sstune),'k:'); hold on;
semilogx(xo,stune - (2*sstune),'k:'); hold on;
semilogx(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('SpatFreq (cyc/deg)');
ylabel('Rate (sp/s)');
set(gca,'XTick',MyTick,'XTickLabel',MyTickLabel);
%*************************************

disp('Computing pref vs non-pref presac');
%***** Re-run analyses, taking only pre-saccadic epoch flashes
%**** run same kind of averages, but broken out by stim type
zkcounts = zeros(1,NTot);
zstcounts = zeros(1,NTot);
zstcounts2 = zeros(1,NTot);
wkcounts = zeros(1,NTot);
wstcounts = zeros(1,NTot);
wstcounts2 = zeros(1,NTot);
vkcounts = zeros(1,NTot);
vstcounts = zeros(1,NTot);
vstcounts2 = zeros(1,NTot);
for k = 1:SN
    if (StimX(k,1))  %only analyze inside range
        stimo = StimX(k,3); % stimulus type
        Sac = 0;
        for j = (k-ENDRNG):(k-STARNG)
           if (j>0) && (j<=SN) 
             if (SacX(j,SO)) 
               Sac = floor( SacX(j,SO) );
               break;
             end
           end
        end
        if (Sac == -1)
         if ( (stimo>0) && (stimo<=NTot))
           if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnt = mean( StimY((k+XX(trepa:trepb))',1)');  % mean over epoch
             zstcounts(stimo) = zstcounts(stimo) + spcnt;
             zstcounts2(stimo) = zstcounts2(stimo) + (spcnt ^ 2);
             zkcounts(stimo) = zkcounts(stimo) + 1;
           end
         end
        end
        if (ismember(Sac,PrefSacs))
         if ( (stimo>0) && (stimo<=NTot))
           if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnt = mean( StimY((k+XX(trepa:trepb))',1)');  % mean over epoch
             wstcounts(stimo) = wstcounts(stimo) + spcnt;
             wstcounts2(stimo) = wstcounts2(stimo) + (spcnt ^ 2);
             wkcounts(stimo) = wkcounts(stimo) + 1;
           end
         end  
        end
        if (ismember(Sac,NonPrefSacs))
          if ( (stimo>0) && (stimo<=NTot))
           if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnt = mean( StimY((k+XX(trepa:trepb))',1)');  % mean over epoch
             vstcounts(stimo) = vstcounts(stimo) + spcnt;
             vstcounts2(stimo) = vstcounts2(stimo) + (spcnt ^ 2);
             vkcounts(stimo) = vkcounts(stimo) + 1;
           end
         end
        end
    end
end
%****
if zkcounts > 0
   zstcounts = zstcounts ./ zkcounts;  % mean
   zstsem = sqrt(  ((zstcounts2 ./ zkcounts) - zstcounts.^2) ./ zkcounts);  % 1 sem 
   zstcounts = zstcounts / SampRate;
   zstsem = zstsem / SampRate;
else
   zstcounts = NaN(size(zstcounts));
   zstsem = NaN(size(zstcounts));
end
%****
if wkcounts > 0
   wstcounts = wstcounts ./ wkcounts;  % mean
   wstsem = sqrt(  ((wstcounts2 ./ wkcounts) - wstcounts.^2) ./ wkcounts);  % 1 sem 
   wstcounts = wstcounts / SampRate;
   wstsem = wstsem / SampRate;
else
   wstcounts = NaN(size(wstcounts));
   wstsem = NaN(size(wstcounts));
end
%****
if vkcounts > 0
   vstcounts = vstcounts ./ vkcounts;  % mean
   vstsem = sqrt(  ((vstcounts2 ./ vkcounts) - vstcounts.^2) ./ vkcounts);  % 1 sem 
   vstcounts = vstcounts / SampRate;
   vstsem = vstsem / SampRate;
else
   zstcounts = NaN(size(zstcounts));
   zstsem = NaN(size(zstcounts));
end

%******* plot out the stimulus selectivity over sig temporal epoch
subplot('position',[0.1 0.05 0.25 0.20]);
uu = reshape(zstcounts,[NOri NSpf]);
su = reshape(zstsem,[NOri NSpf]);
imagesc(SpatOris,SpatFrqs,uu',[mino maxo]); colormap('gray'); h = colorbar;
xlabel('Orientation');
ylabel('Spat Freq');
set(gca,'YTick',MyTick,'YTickLabel',MyTickLabel);
title('Rate(sp/s)');

%***** plot orientation average
%subplot('position',[0.4 0.05 0.225 0.20]);
subplot('position',[0.4 0.35 0.225 0.20]); hold on;
otune = mean(uu');
otune = [otune otune];
sotune = mean(su');
sotune = [sotune sotune];
xo = [SpatOris (SpatOris+180)];
hh = plot(xo,otune,'g-'); hold on;
set(hh,'Linewidth',2);
plot(xo,otune + (2*sotune),'g:'); hold on;
plot(xo,otune - (2*sotune),'g:'); hold on;
plot(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('Orientation (degs)');
ylabel('Rate (sp/s)');

%***** plot spf average
% subplot('position',[0.7 0.05 0.225 0.20]);
subplot('position',[0.7 0.35 0.225 0.20]); hold on;
stune = mean(uu);
sstune = mean(su);
xo = SpatFrqs; 
hh = semilogx(xo,stune,'g-'); hold on;
set(hh,'Linewidth',2);
semilogx(xo,stune + (2*sstune),'g:'); hold on;
semilogx(xo,stune - (2*sstune),'g:'); hold on;
semilogx(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('SpatFreq (cyc/deg)');
ylabel('Rate (sp/s)');
set(gca,'XTick',MyTick,'XTickLabel',MyTickLabel);
%*************************************

%******* plot out the stimulus selectivity over sig temporal epoch
%subplot('position',[0.4 0.05 0.25 0.20]);
uu = reshape(wstcounts,[NOri NSpf]);
su = reshape(wstsem,[NOri NSpf]);
% imagesc(SpatOris,SpatFrqs,uu'); colormap('gray'); h = colorbar;
% xlabel('Orientation');
% ylabel('Spat Freq');
% set(gca,'YTick',MyTick,'YTickLabel',MyTickLabel);
% title('Rate(sp/s)');

%***** plot orientation average
subplot('position',[0.4 0.05 0.225 0.20]);
%subplot('position',[0.4 0.35 0.225 0.20]); hold on;
otune = mean(uu');
otune = [otune otune];
sotune = mean(su');
sotune = [sotune sotune];
xo = [SpatOris (SpatOris+180)];
hh = plot(xo,otune,'r-'); hold on;
set(hh,'Linewidth',2);
plot(xo,otune + (2*sotune),'r:'); hold on;
plot(xo,otune - (2*sotune),'r:'); hold on;
plot(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('Orientation (degs)');
ylabel('Rate (sp/s)');

%***** plot spf average
subplot('position',[0.7 0.05 0.225 0.20]);
%subplot('position',[0.7 0.35 0.225 0.20]); hold on;
stune = mean(uu);
sstune = mean(su);
xo = SpatFrqs; 
hh = semilogx(xo,stune,'r-'); hold on;
set(hh,'Linewidth',2);
semilogx(xo,stune + (2*sstune),'r:'); hold on;
semilogx(xo,stune - (2*sstune),'r:'); hold on;
semilogx(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('SpatFreq (cyc/deg)');
ylabel('Rate (sp/s)');
set(gca,'XTick',MyTick,'XTickLabel',MyTickLabel);
%*************************************

%******* plot out the stimulus selectivity over sig temporal epoch
%subplot('position',[0.7 0.05 0.25 0.20]);
uu = reshape(vstcounts,[NOri NSpf]);
su = reshape(vstsem,[NOri NSpf]);
% imagesc(SpatOris,SpatFrqs,uu'); colormap('gray'); h = colorbar;
% xlabel('Orientation');
% ylabel('Spat Freq');
% set(gca,'YTick',MyTick,'YTickLabel',MyTickLabel);
% title('Rate(sp/s)');

%***** plot orientation average
subplot('position',[0.4 0.05 0.225 0.20]); hold on;
%subplot('position',[0.4 0.35 0.225 0.20]); hold on;
otune = mean(uu');
otune = [otune otune];
sotune = mean(su');
sotune = [sotune sotune];
xo = [SpatOris (SpatOris+180)];
hh = plot(xo,otune,'b-'); hold on;
set(hh,'Linewidth',2);
plot(xo,otune + (2*sotune),'b:'); hold on;
plot(xo,otune - (2*sotune),'b:'); hold on;
plot(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('Orientation (degs)');
ylabel('Rate (sp/s)');

%***** plot spf average
subplot('position',[0.7 0.05 0.225 0.20]); hold on;
% subplot('position',[0.7 0.35 0.225 0.20]); hold on;
stune = mean(uu);
sstune = mean(su);
xo = SpatFrqs; 
hh = semilogx(xo,stune,'b-'); hold on;
set(hh,'Linewidth',2);
semilogx(xo,stune + (2*sstune),'b:'); hold on;
semilogx(xo,stune - (2*sstune),'b:'); hold on;
semilogx(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('SpatFreq (cyc/deg)');
ylabel('Rate (sp/s)');
set(gca,'XTick',MyTick,'XTickLabel',MyTickLabel);
%*************************************

disp('Finished');

return;
