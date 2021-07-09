function PlotFlagGratingKernel(StimX,StimY,FileTag)
% function PlotFlagGratingKernel(Exp)
%   input: takes StimX and StimY (stim history and spike counts)
%     StimX: [Nframes,4] - fields, matlab time, ephys time, stim, (if stim)
%     StimY: [Nframes,1] - spike counts per video frame
%   output: plot of the temporal and tuning kernels
%
% NOTE: VERY SIMILAR TO FORAGE GRATING KERNEL .... build library here?
%

%******** use reverse correlation of stim vs no stim to find temporal
%******** the temporal kernel from the data stream
SampRate = median(diff(StimX(:,1)));  % median screen flip time
TCRIT = 3;  %threshold for significance, in sem
DTA = -5;   % stim frames, relative to stim onset, must be negative
DTB = 25;   % stim frames after onset, must be positive 
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
    % if (StimX(k,4) && StimX(k,5))
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
k = 1;
while k <= size(XX,2)
    if (trepa == 0)  % look for first sig time point
        if ((mcounts(k)-pcounts(k)) > TCRIT*psem(k))
            trepa = k;  % flag first sig time point
        end
    else
        if (trepb == 0)  % find first time point after trepa not sig
            if ((mcounts(k)-pcounts(k)) < TCRIT*psem(k))
                trepb = k;
            end
        end
    end
    k = k + 1;
end
if (trepa > 0) && (trepb == 0)
    trepb = size(XX,2);  % stay sig long term
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
plot(tXX,(pcounts+(TCRIT*psem)),'b:');
plot(tXX,(pcounts-(TCRIT*psem)),'b:');
plot([tXX(1),tXX(end)],[ymean,ymean],'k-');
axis tight;
V = axis;
plot([0,0],[V(3),V(4)],'k-');
xlabel('Time (ms)');
ylabel('Rate (sp/s)');
title(sprintf('Temporal Kernel: %s',FileTag));

%********* if no significant temporal kernel then return ***
if ~trepa || ~trepb
    disp('No significant temporal kernel found, stopping');
    return;
end
plot([tXX(trepa),tXX(trepa)],[V(3),V(4)],'k--');
plot([tXX(trepb),tXX(trepb)],[V(3),V(4)],'k--');

OriSet = unique(StimX(:,3));
if (size(OriSet,1) <= 2)
    disp('Only a single orientation was shown, stopping');
    return;
end

%****** Now using that sig temporal kernel, look for feature selectivity
%****** across spat freq and orientation *********************
SpatOris = (0:11) * 15;   % orientations, figure way to get auto from Exp
NOri = size(SpatOris,2);
SpatFrqs = 4;   % single spat freq, fixed at 4 cycs/deg
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


%***** plot orientation average
subplot('position',[0.1 0.2 0.4 0.3]);
uu = stcounts;
su = stsem;
otune = uu;
otune = [otune otune];
sotune = su;
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

%*************************************

   
return;
