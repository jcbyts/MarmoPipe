function vpx2Ephys = sync_vpx_to_ephys_clock(Exp, ephysTrials, varargin)
% vpx2Ephys = sync_vpx_to_ephys_clock(Exp, ephysTrials)
% Returns a single function handle for converting Eyetracker times into
% the ephys clock times for MarmoV
ip = inputParser();
ip.addParameter('fid', 1)
ip.addParameter('debug', false)
ip.addParameter('mode', 'linear')
ip.parse(varargin{:});

if nargin < 2 || isempty(ephysTrials)
    % Get list of trials with electrophysiolgy timestamps
    ephysTrials = find(cellfun(@(x) ~isnan(x.START_EPHYS), Exp.D));
end

ephysClock = cellfun(@(x) x.START_EPHYS, Exp.D(ephysTrials));
try
    vpxClock = cellfun(@(x) x.START_VPX, Exp.D(ephysTrials));
catch
    vpxClock = cellfun(@(x) x.eyeData(1,6), Exp.D(ephysTrials));
end

ephysClock = [ephysClock; cellfun(@(x) x.END_EPHYS, Exp.D(ephysTrials))];
try
    vpxClock = [vpxClock; cellfun(@(x) x.END_VPX, Exp.D(ephysTrials))];
catch
    vpxClock = [vpxClock; cellfun(@(x) x.eyeData(end,6), Exp.D(ephysTrials))];
end

bad = isnan(ephysClock) | isnan(vpxClock);
ephysClock(bad) = [];
vpxClock(bad) = [];

vpx2Ephys = synchtime.align_clocks(vpxClock, ephysClock, 'mode', ip.Results.mode);

fprintf(1, 'Synchronizing the Ephys and Eye-tracker clocks with %d valid strobes\n', numel(ephysClock));
deltaE = (ephysClock - vpx2Ephys(vpxClock));
totalErrorMs = sum(deltaE.^2)*1e3;

if ip.Results.debug
    figure(999); clf
    subplot(1,2,1)
    plot(vpxClock, ephysClock, '.'); hold on
    plot(vpxClock, vpx2Ephys(vpxClock))

    subplot(1,2,2)
    histogram(deltaE)

    figure(998); clf
    plot(deltaE)
    hold on
    plot(xlim, [0 0], 'k--')
    keyboard
end


assert(totalErrorMs < .1, 'Clock sync failed')
