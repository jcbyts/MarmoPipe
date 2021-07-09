function StimYY = PlotForageGratingKernel(StimX,StimY,FileTag)
% function PlotForageGratingKernel(Exp)
%   input: takes StimX and StimY (stim history and spike counts)
%     StimX: [Nframes,4] - fields, matlab time, ephys time, stim, (if stim)
%     StimY: [Nframes,1] - spike counts per video frame
%   output: plot of the temporal and tuning kernels
%     StimYY: [Nframes,1] - spike counts after removing estimated visual
%                           responses to the oriented grating stimuli
%
%******** use reverse correlation of stim vs no stim to find temporal
%******** the temporal kernel from the data stream
REVTEMPK = false;
%*******
StimYY = StimY;
SampRate = median(diff(StimX(:,1)));  % median screen flip time
DTA = -1;   % stim frames, relative to stim onset, must be negative
DTB = 30;   % stim frames after onset, must be positive 
if (1)
  %*** for MT
  minDTA = ceil(60/(1000 * SampRate))-DTA;  % no earlier than 50 ms
  maxDTA = ceil(120/(1000 * SampRate))-DTA;  % no earlier than 50 ms
else
  %*** for V1
  minDTA = ceil(40/(1000 * SampRate))-DTA;  % no earlier than 50 ms
  maxDTA = ceil(65/(1000 * SampRate))-DTA;
end
%****
fixDTA = [minDTA:maxDTA];
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
   disp('Failed to estimate any temporal kernel, stopping'); 
   return;  % otherwise return, if no temporal kernel there is a problem
end

%******* Find what you consider a significant temporal kernel
trepa = 0;
trepb = 0;
tcrit = 8;  %threshold for significance, in sem
if (~REVTEMPK)
    k = minDTA; %1;
    while k <= size(XX,2)
        if (trepa == 0)  % look for first sig time point
            if (abs(mcounts(k)-pcounts(k)) > tcrit*psem(k)) 
                trepa = k;  % flag first sig time point
            end
            
        else
            if (trepb == 0)  % find first time point after trepa not sig
                if (abs(mcounts(k)-pcounts(k)) < tcrit*psem(k))
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
            if (abs(mcounts(k)-pcounts(k)) > tcrit*psem(k))   
                trepb = k;  % flag first sig time point
            end
        else
            if (trepa == 0)  % find first time point after trepa not sig
                if (abs(mcounts(k)-pcounts(k)) < tcrit*psem(k)) 
                    trepa = k;
                end
            end
        end
        k = k - 1;
    end
end
if (trepa > 0) && (trepb == 0)
    trepb = size(XX,2);  % stay sig long term
end

trepa
trepb

if ~isempty(fixDTA)
    trepa = minDTA;
    trepb = maxDTA;
end

%******* compute summary stats ******
%********* then plot results ********
hf = figure;
set(hf,'position',[100 80 900 900]);
%**** plot the temporal kernel for now
subplot('position',[0.1 0.6 0.8 0.3]);
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
%*********************
tempker = mcounts - ymean;
tempker = tempker / sqrt( sum( tempker .^ 2) );

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
%*********************************************************************
%****** compute the projection of actual temporal kernel on spike counts
%****** for each stim response per ori and spat-freq, gain tuning, and
%****** then you will have a model to de-convolve spike counts
%*******************************************************************
NOri = size(SpatOris,2);
NSpf = size(SpatFrqs,2);
NTot = NOri * NSpf;
%******
gcounts = zeros(1,NTot);
sgcounts = zeros(1,NTot);
sgcounts2 = zeros(1,NTot);
for k = 1:SN
    if (StimX(k,4))  %only analyze inside range
        stimo = StimX(k,3); % stimulus type
        if ( (stimo>0) && (stimo<=NTot))
           if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnt = sum( tempker .* StimY((k+XX),1)' );
             sgcounts(stimo) = sgcounts(stimo) + spcnt;
             sgcounts2(stimo) = sgcounts2(stimo) + (spcnt ^ 2);
             gcounts(stimo) = gcounts(stimo) + 1;
           end
        end
    end
end
if gcounts > 0
   sgcounts = sgcounts ./ gcounts;  % mean
   sgsem = sqrt(  ((sgcounts2 ./ gcounts) - sgcounts.^2) ./ gcounts);  % 1 sem 
else
   sgcounts = NaN(size(sgcounts));
   sgsem = NaN(size(sgcounts));
end
%***********************
hb = figure;
if (0)
  subplot('position',[0.1 0.6 0.4 0.3]);
  plot(tXX,tempker,'k-'); hold on;
  plot(tXX,zeros(size(tempker)),'k--');
end
%subplot('position',[0.55 0.6 0.4 0.4]);
uu = reshape(sgcounts,[NOri NSpf]);
su = reshape(sgsem,[NOri NSpf]);
imagesc(SpatOris,SpatFrqs,uu'); colormap('gray'); h = colorbar;
xlabel('Orientation');
ylabel('Spat Freq');
set(gca,'YTick',MyTick,'YTickLabel',MyTickLabel);
title('Gain projection');

if (0)
  %****** NOW REGRESS OUT THAT PART OF FIRING RATE
  Var_Y = sum( StimY .^ 2);
  StimYY(:,:) = 0; % set to zero first
  for k = 1:SN
    if (StimX(k,4))  %only analyze inside range
        stimo = StimX(k,3); % stimulus type
        if ( (stimo>0) && (stimo<=NTot))
           if ( (k+DTA) > 0) && ((k+DTB)<=SN)
             spcnt = tempker .* sgcounts(stimo); 
             StimYY(k+XX,1) = StimYY(k+XX,1) + spcnt';
           end
        end
    end
  end
  Err_Y = sum( (StimY-StimYY) .^ 2);
  subplot('position',[0.1 0.1 0.8 0.4]);
  plot(1:SN,StimY,'k.-'); hold on;
  plot(1:SN,StimYY,'b-');
  % plot(StimY,StimYY,'k.');
  xlabel('Actual Counts');
  ylabel('Predicted Counts');
  title(sprintf('R2: %6.3f ',(Var_Y - Err_Y)/Var_Y));
  %*********************
end

figure(hf);
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
subplot('position',[0.1 0.1 0.5 0.4]);
uu = reshape(stcounts,[NOri NSpf]);
su = reshape(stsem,[NOri NSpf]);

imagesc(SpatOris,SpatFrqs,uu'); colormap('gray'); h = colorbar;
xlabel('Orientation');
ylabel('Spat Freq');
set(gca,'YTick',MyTick,'YTickLabel',MyTickLabel);
title('Rate(sp/s)');

%***** plot orientation average
subplot('position',[0.7 0.35 0.25 0.15]);
otune = mean(uu');
otune = [otune otune];
sotune = mean(su');
sotune = [sotune sotune];
xo = [SpatOris (SpatOris+180)];
plot(xo,otune,'bo'); hold on;
plot(xo,otune + (2*sotune),'b--'); hold on;
plot(xo,otune - (2*sotune),'b--'); hold on;
plot(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('Orientation (degs)');
ylabel('Rate (sp/s)');

%***** plot spf average
subplot('position',[0.7 0.10 0.25 0.15]);
stune = mean(uu);
sstune = mean(su);
xo = SpatFrqs; 
semilogx(xo,stune,'bo'); hold on;
semilogx(xo,stune + (2*sstune),'b--'); hold on;
semilogx(xo,stune - (2*sstune),'b--'); hold on;
semilogx(xo,ymean*ones(size(xo)),'k-');
axis tight;
V = axis;
mag = 0.1*(V(4)-V(3));
axis([V(1) V(2) (V(3)-mag) (V(4)+mag)]);
xlabel('SpatFreq (cyc/deg)');
ylabel('Rate (sp/s)');
set(gca,'XTick',MyTick,'XTickLabel',MyTickLabel);
%*************************************

   
return;
