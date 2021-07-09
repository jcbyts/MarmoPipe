function [slist] = flag_saccades( eyeSmo , varargin)
% function [slist] = flag_saccades( eyeSmo )
%***
%*** takes a smoothed eye data segment of following format:
%***     size T x 6:   where T is number of samples
%***                   column 1: time (milli-seconds)
%***                   column 2: x position (dva)
%***                   column 3: y position (dva)
%***                   column 4: pupil size
%***                   column 5: velocity of x (dva/s)
%***                   column 6: velocity of y (dva/s)
%*** output:  list of saccades, Kx7 where K is number of saccades
%***   [startsac endsac peaksac kstartsac kendsac kpeaksac flag]];
%***    column 1:  startsac - in time (secs), start of saccade
%***    column 2:  endsac - in time (secs)
%***    column 3:  peaksac - time of peak velocity
%***    column 4:  kstartsac - time as integer of smoEye array
%***    column 5:  kstartsac - time as integer of smoEye array
%***    column 6:  kstartsac - time as integer of smoEye array
%***    column 7:  flag - 0 by default, 4 - if a curved saccades
%***                    and custom analysis can make other distinctions
%***                    like flag the saccade that hits a target
%*******************************************
% Optional Arguments:
% 'VFactor' (default: 5) Relative Velocity Threshold (from Engbert and Trukenbrod, 2014)
% 'MinDuration' (20)     Minumum duration (in ms)
% 'MinGap'      (10)     Minimum gap between saccades before checking if
%                        they are a single saccade (in ms)
% 'FlagCurve'   (1.2)    Flag curved saccades (1 = straight)
% 'SampRate'    (1000)   Sampling rate of the eye trace

% *****************
% Parse optional arguments
ip = inputParser();
ip.KeepUnmatched = true; % don't worry if unmatched arguments are passed in
ip.addParameter('VelThresh', 10)   % dva / s
ip.addParameter('VFactor', 5)      % Relative Velocity Threshold (from Engbert and Trukenbrod, 2014)
ip.addParameter('MinDuration', 20) % in ms
ip.addParameter('MinGap', 10)      % in ms (if gap btw saccades, check if it's only one saccade)
ip.addParameter('FlagCurve', 1.2)  % flag as curved if > this value from straight (1 = straight)
ip.addParameter('SampRate', 1e3)   % sampling rate
ip.parse(varargin{:})

%**************** Parameters for saccade flagging ************
% VelThresh = ip.Results.VelThresh;   % dva/s, threshold on velocity to count as saccade
VFactor   = ip.Results.VFactor;     % velocity factor for flagging saccades
% see Engbert & Trukenbrod, 2014 and their package
MinDuration = ip.Results.MinDuration;  % min duration of saccade (ms)
MinGap      = ip.Results.MinGap;       % If Gap between two flagged saccades is less
% in ms, check if they are really one saccade
FlagCurve   = ip.Results.FlagCurve;    % flag as curved if more than 20% dev from straight
SampRate    = ip.Results.SampRate;
%****************************

Time = eyeSmo(:,1);
Vabs = eyeSmo(:,7);
%***********************
slist = [];

%[VFactor,MinDuration,MinGap,FlagCurve]
%input('in this code');

% main saccade detection routine
[ms,radius] = microsaccMerge(eyeSmo(:,2:3),eyeSmo(:,5:6),VFactor,MinDuration,MinGap);

mylist = [];
if ~isempty(ms)
    ms(:,1:2) = ms(:,1:2) * (1/SampRate);  % covert to secs from ms
    mylist = ms(:,1:2);  % accept flagged saccades as is (Merging already done)
end

%******************************************************************

if ~isempty(mylist)
    %******* fill in full information for slist struct of original program
    for k = 1:size(mylist,1)
        
        %****** compute list of sac properties
        kstartsac = ceil(mylist(k,1) * SampRate);
        kendsac = floor(mylist(k,2) * SampRate);
        startsac = Time(kstartsac);
        endsac = Time(kendsac);
        kband = kstartsac:kendsac;
        %******
        peak = max( Vabs(kband) );
        zz = find( Vabs(kband) == peak );
        kpeaksac = kband(zz(1));
        peaksac = Time(kpeaksac);
        %****
        %[~, zz] = max( Vabs(kband) );
        %kpeaksac = kband(zz(1));
        %peaksac = Time(kpeaksac);
        %****
        %****** finally, check in case it is a curved saccade and mark it
        flag = 0;
        %**** check is first 3/8 move in different vector direction than
        %**** the last 3/8 of the saccade, using start to end saccade time
        kstart = kstartsac;
        kend = kendsac;
        %******* different way to measure curvature:
        %*** 1) measure the distance from saccade start to end
        %*** 2) step along the actual trajectory and integrate its length
        %*** 3) then ask if actual path is within 25% of the straight path
        svec = [ (eyeSmo(kend,2) - eyeSmo(kstart,2)),...
            (eyeSmo(kend,3) - eyeSmo(kstart,3))];
        straightpath = norm(svec);
        integratepath = 0;
        for zk = kstart:(kend-1)
            svec2 = [ (eyeSmo(zk+1,2) - eyeSmo(zk,2)),...
                (eyeSmo(zk+1,3) - eyeSmo(zk,3))];
            integratepath = integratepath + norm(svec2);
        end
        if (integratepath > (FlagCurve * straightpath) )
            flag = 4;
        end
        %*****************************
        slist = [slist ; [startsac endsac peaksac kstartsac kendsac kpeaksac flag]];
    end
end

