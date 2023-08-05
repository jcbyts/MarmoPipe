function ptb2Ephys = sync_ptb_to_ephys_clock(Exp, ephysTrials, varargin)
% ptb2Ephys = syncPtb2EphysClock(Exp)
% Returns a single function handle for converting Psychtoolbox times into
% the ephys clock times for MarmoV
ip = inputParser();
ip.addParameter('fid', 1)
ip.addParameter('debug', false)
ip.addParameter('mode', 'linear')
ip.parse(varargin{:});

fid = ip.Results.fid;

if nargin < 2 || isempty(ephysTrials)
    % Get list of trials with electrophysiolgy timestamps
    ephysTrials = find(cellfun(@(x) ~isnan(x.START_EPHYS), Exp.D));
end

ephysClock = cellfun(@(x) x.START_EPHYS, Exp.D(ephysTrials));
try
    ptbClock = cellfun(@(x) x.STARTCLOCKTIME, Exp.D(ephysTrials));
catch
    ptbClock = cellfun(@(x) x.eyeData(1,6), Exp.D(ephysTrials));
end

ephysClock = [ephysClock; cellfun(@(x) x.END_EPHYS, Exp.D(ephysTrials))];
try
    ptbClock = [ptbClock; cellfun(@(x) x.ENDCLOCKTIME, Exp.D(ephysTrials))];
catch
    ptbClock = [ptbClock; cellfun(@(x) x.eyeData(end,6), Exp.D(ephysTrials))];
end


bad = isnan(ephysClock) | isnan(ptbClock);
ephysClock(bad) = [];
ptbClock(bad) = [];

ptb2Ephys = synchtime.align_clocks(ptbClock, ephysClock, 'mode', ip.Results.mode);

fprintf(fid, 'Synchronizing the Ephys and PTB clocks with %d valid strobes\n', numel(ephysClock));
deltaE = (ephysClock - ptb2Ephys(ptbClock));
totalErrorMs = sum(deltaE.^2)*1e3;
fprintf(fid, 'Total error (SSE): %02.5f ms\n', totalErrorMs);

if ip.Results.debug
    figure(999); clf
    subplot(1,2,1)
    plot(ptbClock, ephysClock, '.'); hold on
    plot(ptbClock, ptb2Ephys(ptbClock))

    subplot(1,2,2)
    histogram(deltaE)

    figure(998); clf
    plot(deltaE)
    hold on
    plot(xlim, [0 0], 'k--')
    keyboard
end


assert(totalErrorMs < .1, 'Clock sync failed')
