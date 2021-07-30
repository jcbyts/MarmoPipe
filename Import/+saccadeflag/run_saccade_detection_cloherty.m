function Exp = run_saccade_detection_cloherty(Exp, varargin)
% RUN SACCADE DETECTION
% Detects saccades for each trial of an Exp struct using the algorithm
% described in Cloherty,Yates et al. (2020)
%
%-------------------------------------------------------------------
% Input:
%
% Exp [struct]: the marmoview struct
%               must contain fields vpx and D
%
% Output:
%
% Exp [struct]: modified marmoview struct contains slist
%   slist [array]: K x 7, where K is number of saccades
%       [startsac endsac peaksac kstartsac kendsac kpeaksac flag]]
%           column 1:  startsac - in time (secs), start of saccade
%           column 2:  endsac - in time (secs)
%           column 3:  peaksac - time of peak velocity
%           column 4:  kstartsac - time as integer of smoEye array
%           column 5:  kstartsac - time as integer of smoEye array
%           column 6:  kstartsac - time as integer of smoEye array
%           column 7:  flag - 0 by default, 4 - if a curved saccades
%                      and custom analysis can make other distinctions
%                      like flag the saccade that hits a target
%
% Optional Arguments:
%
% 'VFactor'     (5)     Relative Velocity Threshold (from Engbert and Trukenbrod, 2014)
% 'MinDuration' (20)    Minumum duration (in ms)
% 'MinGap'      (10)    Minimum gap between saccades before checking if
%                       they are a single saccade (in ms)
% 'FlagCurve'   (1.2)   Flag curved saccades (1 = straight)
% 'SampRate'    (1000)  Sampling rate of the eye trace

% see also: microsaccMerge, +saccadeflag.flag_saccades, +saccadeflag.eyeposition_smooth


ip = inputParser();
ip.KeepUnmatched = true;
ip.addParameter('accthresh', 2e4)
ip.addParameter('velthresh', 8)
ip.addParameter('velpeak', 10)
ip.addParameter('isi', 0.04)
ip.addParameter('ShowTrials', false)
ip.addParameter('fid', 1)
ip.parse(varargin{:});

% % convert unmatched arguments for input to the saccade_flag function
% args = [fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]';
fid = ip.Results.fid;
ShowTrials = ip.Results.ShowTrials;

if ShowTrials
    H = figure(1); clf
    set(H,'position',[100 200 800 800]);
end


vtt = Exp.vpx.smo(:,1);
vxxd = Exp.vpx.smo(:,2);
vyyd = Exp.vpx.smo(:,3);
vxx = medfilt1(vxxd, 5);
vyy = medfilt1(vyyd, 5);
vxx = imgaussfilt(vxx, 7);
vyy = imgaussfilt(vyy, 7);

% vpp = Exp.vpx.smo(:,4);
% vx = Exp.vpx.smo(:,5);
% vy = Exp.vpx.smo(:,6);
amp = hypot(Exp.vpx.smo(:,2),Exp.vpx.smo(:,3));
spd = [0; diff(amp)];
Exp.vpx.smo(:,7) = spd;
% spd = ;



goodix = find(amp < 20);
saccades = saccadeflag.find_saccades(Exp.vpx.smo(goodix,1),vxx(goodix),vyy(goodix),...
    'accthresh', ip.Results.accthresh,...
    'velthresh', ip.Results.velthresh,...
    'velpeak', ip.Results.velpeak,...
    'isi', ip.Results.isi, ...
    'debug', false);

% handle mouse input
if isempty(saccades)
    slist = [];
else
    slist = [saccades.tstart, saccades.tend, saccades.tpeak, ...
    goodix(saccades.startIndex), goodix(saccades.endIndex), goodix(saccades.peakIndex), zeros(numel(saccades.tstart), 1)];

    % get proper times for the saccades
    slist(:,1) = vtt(slist(:,4));
    slist(:,2) = vtt(slist(:,5));
    slist(:,3) = vtt(slist(:,6));
end

% store values
Exp.slist = slist;


% --- Saccade / Eye Trace QA
nTimePoints = numel(vxxd);
Exp.vpx.Labels = zeros(nTimePoints, 1);
Exp.vpx.Labels(:) = 1; % initialize everything to fixation
Exp.vpx.LabelIds = {'Fixation', 'Saccade', 'Blink', 'Lost'};

% insert saccades
nSaccades = size(Exp.slist,1);
fprintf(fid, 'Found %d saccades\n', nSaccades);

if nSaccades == 0
    return
end
% 
% expand blinks with some padding

blink = Exp.vpx.smo(:,4)==0; % when pupil size is zero
if all(blink) % means eye position was acquired online. this is meaningless
    blink(:) = false;
end
padding = 5;
bc = ones(padding,1);
blink = filtfilt(bc, 1, double(blink)) > 0;
fprintf(fid, '%d/%d samples excluded because pupil size was 0\n', sum(blink), numel(blink));
Exp.vpx.Labels(blink) = 4;

% flag velocities that make no sense
badVel = Exp.vpx.smo(:,7) > 1e3; % unphysiologically plausible  velocities 
Exp.vpx.Labels(badVel) = 4;
v = Exp.vpx.smo(Exp.slist(:,6),7);
fprintf(fid, '%d/%d samples excluded because eye jumped > 1000 deg/sec\n', sum(v > 1e3), numel(v));
Exp.slist(v > 1e3,:) = [];

% tracker is unreliable more than 20 d.v.a from center (flag and label)
offScreen = hypot(Exp.vpx.smo(:,2),Exp.vpx.smo(:,2)) > 20;
fprintf(fid, '%d/%d samples excluded because eye pos was off screen (> 20 deg)\n', sum(offScreen), numel(offScreen));
Exp.vpx.Labels(offScreen) = 4;

% if any saccades intersect with blinks, remove them
bad = blink(Exp.slist(:,4)) | blink(Exp.slist(:,5));
Exp.slist(bad,:) = [];

if isfield(Exp.vpx, 'raw')
    outlier = isnan(Exp.vpx.raw(:,2));
    Exp.vpx.Labels(outlier) = 4;
    bad = outlier(Exp.slist(:,4)) | outlier(Exp.slist(:,5));
    Exp.slist(bad,:) = [];
end

% insert saccades
nSaccades = size(Exp.slist,1);
fprintf(fid, 'Ended with %d saccades after QA\n', nSaccades);
for iSac = 1:nSaccades
    ix = Exp.slist(iSac,4):Exp.slist(iSac,5);
    Exp.vpx.Labels(ix) = 2;
end

%% plot it
if ShowTrials
    
   Labels = Exp.vpx.Labels;
    slist = Exp.slist;
    for iSac = 1:nSaccades
        ix = slist(iSac,4):slist(iSac,5);
        Labels(ix) = 2;
    end
    
    %%
    figure(1); clf
    for  i = 1:2
        subplot(2,1,1)
        ix = find(Labels==i);
        plot(ix, vxx(ix), '.'); hold on
        plot(ix, vyy(ix), '.');
        
        subplot(2,1,2)
        plot(ix, spd(ix), '.'); hold on
        
    end
    drawnow
    
    
    
    
    figure(2); clf
    xs = Exp.vpx.smo(Exp.slist(:,4),2);
    ys = Exp.vpx.smo(Exp.slist(:,4),3);
    xe = Exp.vpx.smo(Exp.slist(:,5),2);
    ye = Exp.vpx.smo(Exp.slist(:,5),3);

    dx = xe-xs;
    dy = ye-ys;
    dd =  hypot(dx, dy);

    dt = Exp.vpx.smo(Exp.slist(:,5),1) - Exp.vpx.smo(Exp.slist(:,4),1);
    v = Exp.vpx.smo(Exp.slist(:,6),7);
    
    subplot(1,2,1)
    plot(dd, v, '.')
    
    subplot(1,2,2)
    plot(dd, dt, '.')
    
    win = [-100 500];
    st = Exp.slist(:,4);
    en = Exp.slist(:,5);
    st = st(2:end);
    en = en(1:end-1);
    dur = st - en;
    [~, ind] = sort(dur);
    [an, ~, ~, wfs] = eventTriggeredAverage(spd, en(ind), win);
    figure, imagesc(log(wfs))
    
%     figure, imagesc(log(wfs(1:100,:)))
end

