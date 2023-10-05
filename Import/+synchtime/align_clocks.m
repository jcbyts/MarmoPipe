function [syncFun, invSyncFun] = align_clocks(clock1, clock2, varargin)
% syncFun = align_clocks(clock1, clock2, 'mode', MODE)

ip = inputParser();
ip.addParameter('mode', 'linear')
ip.parse(varargin{:});

bad = isnan(clock1) | isnan(clock2);
clock1(bad) = [];
clock2(bad) = [];

switch ip.Results.mode
    
    case 'lsq'
        % % least-squares to synchronize
        X = [clock2 ones(numel(clock1), 1)];
        w = (X'*X)\(X'*clock1);
        syncFun = @(t) (t - w(1))/w(2);
        invSyncFun = @(t) t*w(2) + w(1);

    case 'robust'
        % robust least-squares to synchronize
        w = robustfit(clock2, clock1);
        syncFun = @(t) (t - w(1))/w(2);
        invSyncFun = @(t) t*w(2) + w(1);
        

    case {'linear', 'pchip', 'spline', 'cubic', 'makima'}
        syncFun = @(t) interp1(clock1, clock2, t, ip.Results.mode, 'extrap');
        invSyncFun = @(t) interp1(clock2, clock1, t, ip.Results.mode, 'extrap');
    
    otherwise
        error('align_clocks: unrecognized method for aligning clocks')
end

fprintf('Synchronizing Clock 1 and Clock 2 with %d timestamps using [%s]\n', numel(clock2), ip.Results.mode)
pred = syncFun(clock1);
totalErrorMs = sum((clock2 - pred).^2)*1e3;
fprintf('Total error (SSE): %02.5f ms\n', totalErrorMs)