function fit = fit_vonmises(theta, count, useBootstrapping)
% fit = fit_vonmises(th, R)
% fit vonmises with proper min and max using maximum likelihood with
% poisson distributed spike counts
% Assumes directions are in degrees
if nargin < 3
    useBootstrapping = false;
end

% initialize parameters
thetas = unique(theta);

tuningCurve = arrayfun(@(x) mean(count(theta==x)), thetas);
tuningCurveSE = arrayfun(@(x) std(count(theta==x))./sqrt(sum(theta==x)), thetas);

% initialize parameters
minFR = min(tuningCurve);
maxFR = max(tuningCurve);
r = sum(count.*exp(1i*(theta/180*pi)))/sum(count); %resultant length
thPref = wrapTo2Pi(angle(r)); % circular mean
Kappa  = (1 - abs(r)); % circular variance

params0 = [minFR maxFR, Kappa, thPref];

opts = optimset('MaxFunEval', 10e3, 'Display', 'off');

% if you want to try bounds
LB = [0 0 0.01 0]; % lower bound
UB = [maxFR maxFR*1.5 50 2*pi];

% LB = [];
% UB = [];
if useBootstrapping
    nBoot = 100; %2000;
    pBoot = zeros(nBoot,numel(params0));
    
    nTrials = numel(count);
    
    for i = 1:nBoot
        inds = randi(nTrials, nTrials, 1);
        fun = @(params) vonMises.neglogli_poissGLM(vonMises.vonmises(theta(inds)/180*pi, params), count(inds));
        pBoot(i,:) = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);
    end
    
    paramsSD = prctile(pBoot, [16 84]); % 68% confidence intervals
    
    %build objective function
    fun = @(params) vonMises.neglogli_poissGLM(vonMises.vonmises(theta/180*pi, params), count);
    
    % optimization
    [phat, fval, ~, ~, ~, ~, H] = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);
else
    %build objective function
    fun = @(params) vonMises.neglogli_poissGLM(vonMises.vonmises(theta/180*pi, params), count);
    
    % optimization
    [phat, fval, ~, ~, ~, ~, H] = fmincon(fun, params0, [], [], [], [], LB, UB, [], opts);
    
    % error bars
    paramsSD = sqrt(diag(inv(H)))';
end

% fit output
fit.paramsML = phat;
fit.fvalue = fval;
fit.paramsSD = paramsSD;
fit.tuningFun = @(th) vonMises.vonmises(th/180*pi, phat);
fit.vonmises = @(th, params) vonMises.vonmises(th/180*pi, params);
if exist('pBoot', 'var')
    fit.pBoot = pBoot;
    
    thetaPref = pBoot(:,4)/pi*180;
    halfWidth = k2hw(pBoot(:,3))/pi*180;
    
    fit.thetaPref = phat(4)/pi*180;
    fit.halfWidth = k2hw(phat(3))/pi*180;
    
    cibnd = [16 84];
    fit.thetaPrefSD = abs(prctile(thetaPref, cibnd) - fit.thetaPref);
    fit.halfWidthSD = abs(prctile(halfWidth, cibnd) - fit.halfWidth);
    fit.minFR = phat(1);
    fit.minFRSD = abs(prctile(pBoot(:,1), cibnd) - phat(1));
    fit.maxFR = phat(2);
    fit.maxFRSD = abs(prctile(pBoot(:,2), cibnd) - phat(2));
else
    fit.thetaPref = phat(4)/pi*180;
    fit.thetaPrefSD = paramsSD(4)/pi*180;
    fit.halfWidth = k2hw(phat(3))/pi*180;
    fit.halfWidthSD = abs(k2hw(phat(3)+paramsSD(3))/pi*180 - fit.halfWidth);
    fit.minFR = phat(1);
    fit.minFRSD = paramsSD(1);
    fit.maxFR = phat(2);
    fit.maxFRSD = paramsSD(2);
end

fit.thetas = thetas;
fit.tuningCurve = tuningCurve;
fit.tuningCurveSE = tuningCurveSE;

function hw = k2hw(k)
% von Mises k to half-width at half-maximum (in rad.)
% hw = acos(log(0.5)./k+1);
hw = acos(log(.5 + .5*exp(2*k))./k -1);

