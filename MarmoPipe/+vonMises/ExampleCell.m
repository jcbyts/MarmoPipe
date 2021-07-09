%% load data
load SampleData.mat

cond = SampleData(:,1);
theta = SampleData(:,2);
count = SampleData(:,3);

useBootstrapping = true;
% fit using bootstrapped errorbars
fit1 = fit_vonmises(theta(cond==1), count(cond==1), useBootstrapping);
fit2 = fit_vonmises(theta(cond==0), count(cond==0), useBootstrapping);

%% plotting

if useBootstrapping
    
    th = linspace(0, 360, 100); % directions for plotting the fits
    
    % generate tuning curves from each bootstrap parameter set
    nBoot = size(fit1.pBoot,1);
    tc1 = nan(nBoot,numel(th));
    tc2 = nan(nBoot,numel(th));
    for i = 1:nBoot
        tc1(i,:) = fit1.vonmises(th, fit1.pBoot(i,:));
        tc2(i,:) = fit2.vonmises(th, fit2.pBoot(i,:));
    end
    
    % use the percentiles of the output curves to generate confidence
    % intervals on the tuning curves
    
    % ci = [2.5 97.5];
    ci = [16 84]; % 68% cofidence intervals 
    figure(1); clf
    %     subplot(1,2,1)
    
    ix = cond==1;
    plot(theta(ix), count(ix), '.r'); hold on
    plot(th, mean(tc1), 'r'); hold on
    plot(th, prctile(tc1, ci), 'r--')
    
    ix = cond==0;
    plot(theta(ix)+.1, count(ix), '.b');
    plot(th, mean(tc2), 'b');
    plot(th, prctile(tc2, ci), 'b--')
    set(gca, 'XTick', 0:45:360, 'TickDir', 'out', 'box', 'off')
    xlabel('Direction')
    ylabel('Spike Count')
    
    % plot parameters
    %     subplot(1,2,2)
    figure(2); clf
    offset = .1;
    errorbar(1, fit1.thetaPref, fit1.thetaPrefSD(1),fit1.thetaPrefSD(2), 'or', 'MarkerFaceColor', 'r', 'Capsize', 0); hold on
    errorbar(1+offset, fit2.thetaPref, fit2.thetaPrefSD(1),fit2.thetaPrefSD(2), 'ob', 'MarkerFaceColor', 'b', 'Capsize', 0); hold on
    
    errorbar(2, fit1.halfWidth, fit1.halfWidthSD(1),fit1.halfWidthSD(2), 'or', 'MarkerFaceColor', 'r', 'Capsize', 0); hold on
    errorbar(2+offset, fit2.halfWidth, fit2.halfWidthSD(1),fit2.halfWidthSD(2), 'ob', 'MarkerFaceColor', 'b', 'Capsize', 0); hold on
    
    errorbar(3, fit1.minFR, fit1.minFRSD(1),fit1.minFRSD(2), 'or', 'MarkerFaceColor', 'r', 'Capsize', 0); hold on
    errorbar(3+offset, fit2.minFR, fit2.minFRSD(1),fit2.minFRSD(2), 'ob', 'MarkerFaceColor', 'b', 'Capsize', 0); hold on
    
    errorbar(4, fit1.maxFR, fit1.maxFRSD(1),fit1.maxFRSD(2), 'or', 'MarkerFaceColor', 'r', 'Capsize', 0); hold on
    errorbar(4+offset, fit2.maxFR, fit2.maxFRSD(1),fit2.maxFRSD(2), 'ob', 'MarkerFaceColor', 'b', 'Capsize', 0); hold on
    
    
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Pref. Dir.', 'HWHM', 'MinFR', 'MaxFR'}, 'XTickLabelRotation', -45)
    
    set(1, 'PaperSize', [4 4], 'PaperPosition', [0 0 4 4])
    
else
    
    figure(1); clf
    subplot(1,2,1)
    
    th = linspace(0, 360, 100); % directions for plotting the fits
    
    % towards condition
    ix = cond==1;
    plot(theta(ix)+.1, count(ix), '.r'); hold on
    plot(th, fit1.tuningFun(th), 'r'); hold on
    
    warning('this does not work for proper tuning curve errorbars...')
    plot(th, fit1.vonmises(th, fit1.paramsML + fit1.paramsSD), 'r--')
    plot(th, fit1.vonmises(th, fit1.paramsML - fit1.paramsSD), 'r--')
    
    % away condition
    ix = cond==0;
    plot(theta(ix), count(ix), '.b');
    plot(th, fit2.tuningFun(th), 'b')
    plot(th, fit2.vonmises(th, fit2.paramsML + fit2.paramsSD), 'b--')
    plot(th, fit2.vonmises(th, fit2.paramsML - fit2.paramsSD), 'b--')
    
    subplot(1,2,2)
    errorbar(fit1.paramsML([4 3 1 2]), fit1.paramsSD([4 3 1 2]), 'or'); hold on
    errorbar(fit2.paramsML([4 3 1 2]), fit2.paramsSD([4 3 1 2]), 'ob'); hold on

    set(gca, 'XTick', 1:4, 'XTickLabel', {'Pref. Dir. (rad)', 'kappa', 'MinFR', 'MaxFR'}, 'XTickLabelRotation', -45)
end