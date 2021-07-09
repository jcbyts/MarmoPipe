function Exp = run_saccade_detection_trial(Exp, varargin)
% RUN SACCADE DETECTION
% Detects saccades for each trial of an Exp struct using the algorithm
% described in Engbert and Mergenthaler (2006)
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
ip.addParameter('ShowTrials', false)
ip.parse(varargin{:});

% convert unmatched arguments for input to the saccade_flag function
args = [fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]';

ShowTrials = ip.Results.ShowTrials;

if ShowTrials
    H = figure(1); clf
    set(H,'position',[100 200 800 800]);
end


%***** loop over trials
for iTrial = 1:size(Exp.D,1)
    
    %******* check for timing strobe
    sVPX = Exp.D{iTrial}.START_VPX;
    eVPX = Exp.D{iTrial}.END_VPX;
    
    if (~isnan(sVPX) && ~isnan(eVPX) )  % make sure it was recorded OK
        
        %********** matlab eye position (online)
        tt = Exp.D{iTrial}.eyeData(:,1) - Exp.D{iTrial}.eyeData(6,1); % sub first time
        xx = Exp.D{iTrial}.eyeData(:,2);
        yy = Exp.D{iTrial}.eyeData(:,3);
        
        %********* get VPX eye position
        zz = find( (Exp.vpx.smo(:,1) >= sVPX) & (Exp.vpx.smo(:,1) < eVPX) );
        if isempty(zz)
            fprintf('Skipping trial %d\n', iTrial)
            continue
        end
        vtt = Exp.vpx.smo(zz,1) - Exp.vpx.smo(zz(1),1);  %subtract first time
        vxx = Exp.vpx.smo(zz,2);
        vyy = Exp.vpx.smo(zz,3);
        vpp = Exp.vpx.smo(zz,4);
        
        %**** eye transform
        pixPerDeg = Exp.S.pixPerDeg;
        c = Exp.D{iTrial}.c;
        dx = Exp.D{iTrial}.dx;
        dy = Exp.D{iTrial}.dy;
        
        %******transform positions
        xx = (xx - c(1))/(dx*pixPerDeg);
        yy = (yy - c(2))/(dy*pixPerDeg);
        
        %*** same on VPX data
        vxx = (vxx - c(1))/(dx*pixPerDeg);
        vyy = 1 - vyy;
        vyy = (vyy - c(2))/(dy*pixPerDeg);
        
        %******** process eye position to smooth and comp vel
        eye = [vtt,vxx,vyy,vpp];
        eyeSmo = +saccadeflag.eyeposition_smooth(eye);
        Exp.D{iTrial}.eyeSmo = eyeSmo;
        
        if any(isnan(eyeSmo(1,:)))
            eyeSmo(1,:) = zeros(1,7);
        end
        
        for i = 1:size(eyeSmo,2)
            try
            eyeSmo(:,i) = repnan(eyeSmo(:,i), 'nearest');
            end
        end
        %******** process smoothed eye to flag saccades
        slist = +saccadeflag.flag_saccades(eyeSmo, args{:});
        
        %********* store new saccade list
        Exp.D{iTrial}.slist = slist;
        
        %*********************************************
        
        %********* plot eye pos over time
        if ShowTrials
            subplot('position',[0.075 0.575 0.85 0.35]); hold off;
            plot(tt,xx,'r:'); hold on;
            plot(tt,yy,'b:'); hold on;
            legend('H','V');
            plot(vtt,vxx,'r.'); hold on;
            plot(vtt,vyy,'b.'); hold on;
            %**** plot smoothed points
            plot(eyeSmo(:,1),eyeSmo(:,2),'r-');
            plot(eyeSmo(:,1),eyeSmo(:,3),'b-');
            V = axis;
            if (~isempty(slist))
               for k = 1:size(slist,1)
                   if (slist(k,7) == 1)
                       colo = 'g';
                   else
                       if (slist(k,7) == 4)
                           colo = 'm';
                       else
                           colo = 'k';
                       end
                   end
                   plot([slist(k,1),slist(k,1)],[V(3),V(4)],[colo,'-']);
                   plot([slist(k,2),slist(k,2)],[V(3),V(4)],[colo,'-']);
               end
            end
            %*******
            title(sprintf('Trial = %d',iTrial));
            xlabel('Time (ms)');
            ylabel('EyePos');
            %*******
            if isfield(Exp.D{iTrial}.PR,'error')
                Error = Exp.D{iTrial}.PR.error;
            else
                Error = 0;
            end
            %****************
            subplot('position',[0.1 0.075 0.4 0.4]); hold off;
            plot(vxx,vyy,'k.'); hold on;
            plot(eyeSmo(:,2),eyeSmo(:,3),'k-');
            plot(0,0,'k+');
            axis([-10 10 -10 10]);
            title(sprintf('Error = %d',Error));
            xlabel('X Pos)');
            ylabel('Y Pos');    
        end 
        %****************
        if ShowTrials
           input('check');
        else
           fprintf('Processing trial %d\n',iTrial);
        end
        %************
    else
        Exp.D{iTrial}.eyeSmo = [];  % if not recorded, put empty record
        Exp.D{iTrial}.slist = []; 
    end
end  % loop over all trials