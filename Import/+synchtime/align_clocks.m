function [syncFun, invSyncFun] = align_clocks(clock1, clock2)
% syncFun = align_clocks(clock1, clock2)

bad = isnan(clock1) | isnan(clock2);
clock1(bad) = [];
clock2(bad) = [];

w = robustfit(clock2, clock1);
% 
% % least-squares to synchronize
% X = [clock2 ones(numel(clock1), 1)];
% w = (X'*X)\(X'*clock1);
% syncFun = @(t) (t - w(2))/w(1);

% function to synchronize
syncFun = @(t) (t - w(1))/w(2);
invSyncFun = @(t) t*w(2) + w(1);

fprintf('Synchronizing Clock 1 and Clock 2 with %d timestamps\n', numel(clock2))
pred = syncFun(clock1);
totalErrorMs = sum((clock2 - pred).^2)*1e3;
fprintf('Total error (SSE): %02.5f ms\n', totalErrorMs)

% remove outliers, re-run
adiff = abs(clock2 - pred);
iix = (adiff / median(adiff)) < 20;

totalErrorMs = sum((clock2(iix) - syncFun(clock1(iix))).^2)*1e3;
fprintf('Total error (SSE) without outliers: %02.5f ms\n', totalErrorMs)
fprintf('ignoring %d/%d points\n', sum(~iix), numel(iix))