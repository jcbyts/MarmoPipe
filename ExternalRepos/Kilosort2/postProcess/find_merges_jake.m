function rez = find_merges(rez, flag)
% This is not the publically available Kilosort merge script

% Jake wrote this because it roughly matches the heuristics he follows
% while sorting.


% 1. Triage
%       The first step is to triage potential single units from obvious
%       garbage. One reasonable way to rule out SU is if the cluster has
%       tons of refactory period violations. We can estimate the cluster
%       firing rate (if it were a homogenous poisson process) and then
%       check the rate within the refractory period. That ratio
%       (fpRate/totalRate) should be really close to 0 if the cell is a SU.
% 
ops = rez.ops;
dt = 1/1000;

Xsim = rez.simScore;
Nk = size(Xsim,1);
Xsim = Xsim - diag(diag(Xsim));

tfi = rez.iNeigh;
tf = rez.cProj;

refDur = ceil(rez.ops.fs/1e3); % 1ms in samples

temps = reshape(rez.Wraw, [], Nk);
wfCorr = corr(temps); % how correlated are the unwhitened templates
wfAmp = squeeze(mean(abs(rez.Wraw(:,:,:)),2));
wfAmp = wfAmp.^2 ./ sum(wfAmp.^2);
wfCtr = rez.yc'*wfAmp;


% sort by firing rate first
nspk = zeros(Nk, 1);
fpRate = zeros(Nk,1);
for j = 1:Nk
    ix = rez.st3(:,2)==j; % index for this unit
    nspk(j) = sum(ix);
    st = sort(rez.st3(ix,1));
    T = st(end)-st(1);
    totalRate = nspk(j)/T;
    nrefviol = sum(abs(diff(st)) < refDur);
    refRate = nrefviol / (refDur * nspk(j));
    fpRate(j) = min( refRate / totalRate, 1);
end

putativeSUs = find(fpRate < .2);

[~, isort] = sort(nspk);
fprintf('initialized spike counts\n')

if ~flag
   rez.R_CCG = Inf * ones(Nk);
   rez.Q_CCG = Inf * ones(Nk);
   rez.K_CCG = {};
end


for j = 1:Nk
    fprintf('cluster %d\n', j)
    jIsSU = ismember(isort(j), putativeSUs);
    s1 = rez.st3(rez.st3(:,2)==isort(j), 1)/ops.fs;
    if numel(s1)~=nspk(isort(j))
        fprintf('lost track of spike counts')
    end    
    
    % select possible merges
    % forget the similarity matrix. Only consider neurons on the same group
    % of channels and with the same SU status. This is the main triage
    distance = sqrt( (wfCtr(isort(j)) - wfCtr).^2);
    distThresh = 400;
    corrThresh = 0.4;
    
    possible = find((distance < distThresh) & (wfCorr(isort(j),:) > corrThresh));
    if jIsSU
        possible = intersect(possible, putativeSUs);
    else
        possible = intersect(possible, setdiff(1:Nk, putativeSUs));
    end
    possible = setdiff(possible, isort(j)); % obviously remove self
    
    [~, xind] = sort(Xsim(isort(j),possible), 'descend');
    possible = possible(xind); % units we will consider for merge, sorted by their similarity score
    
%     [ccsort, ix] = sort(Xsim(isort(j),:) .* (nspk'>numel(s1)), 'descend');
%     ienu = find(ccsort<0.1,1) - 1;    
    for k = possible(:)'
        
        kIsSU = ismember(k, putativeSUs);
        assert( kIsSU == jIsSU, 'find_merges: something wrong. k and j must be of same SU status')
        
        ix1 = rez.st3(:,2)==isort(j);
        ix2 = rez.st3(:,2)==k;
        
        if 1
            f1 = mean(tf(ix1,:));
            f1 = f1(:) / norm(f1);
            
            f2 = mean(tf(ix2,:));
            f2 = f2(:) / norm(f2);
            
            jf1 = tf(ix1,:)*f1;
            jf2 = tf(ix1,:)*f2;
            kf1 = tf(ix2,:)*f1;
            kf2 = tf(ix2,:)*f2;
            
            
            m = mean([jf1 jf2]) - mean([kf1 kf2]);
            X = [[jf1 jf2]; [kf1 kf2]];
            C = cov(X,1);
            w = C\m';
            
            % project on weights for visualization
            x1 = [jf1 jf2]*w;
            x2 = [kf1 kf2]*w;
            
            % percent correct discrimination under gaussian approximation
            p = normcdf(sqrt(m*w)/2);
            
            % plot outcome
            figure(1); clf
            subplot(1,2,1) % waveform templates for both units
            yoff = (1:numel(rez.yc))';
            plot(bsxfun(@plus, rez.Wraw(:,:,isort(j)), yoff)', 'b'); hold on
            plot(bsxfun(@plus, rez.Wraw(:,:,k), yoff)', 'r');
            title(wfCorr(isort(j),k))
            
            subplot(3,2,2) % projection on PCs
            plot(kf1, kf2, 'r.'); hold on
            plot(jf1, jf2, 'b.');
            
            subplot(3,2,4) % discriminability of PC space
            h = histogram(x2, 'Normalization', 'probability'); hold on
            histogram(x1, 'BinEdges', h.BinEdges, 'Normalization', 'probability');
            
            title(p)
            
            s2 = rez.st3(rez.st3(:,2)==k, 1)/ops.fs;
            [K, Qi, Q00, Q01, rir] = ccg(s1, s2, 500, dt);
            
            subplot(3,2,6) % CCG
            plot(K)
            Q = min(Qi/(max(Q00, Q01)));
            
            R = min(rir);
            title([Q R])
            drawnow
            doMerge = (p < .9) | (Q < .3 & R < .5);
        end
            
        
%         pause
        if flag
            if doMerge
%             if Q<.2 % && R<.05
                
                % do not merge simultaneous spikes, it doesn't make sense
                % if both clusters spiked at the same time, then that should
                % become one spike, not two
                ti = rez.st3(rez.st3(:,2)==k,1);
                tjinds = find(rez.st3(:,2)==isort(j));
                tj = rez.st3(tjinds,1);
                nj = numel(tj);
                viol = false(nj,1);
                for ispk = 1:nj
                   viol(ispk) = any(abs(tj(ispk) - ti) < refDur);
                end
                
                % now merge j into i and move on
                ijmerge = rez.st3(:,2)==isort(j); % find spike times to merge
                ijmerge(tjinds(viol)) = false; % ignore simultaneous spikes
                rez.st3(tjinds(viol),2) = Nk + 1; % make the violations an invalid cluster
                rez.st3(ijmerge,2) = k;
                
                nspk(k) = nspk(k) + nspk(isort(j));
                fprintf('merged %d into %d \n', isort(j), k)
                break; % break loop because we've now merged j into k
            end
        else
            rez.R_CCG(isort(j), ix(k)) = R;
            rez.Q_CCG(isort(j), ix(k)) = Q;
            
            rez.K_CCG{isort(j), ix(k)} = K;                        
            rez.K_CCG{ix(k), isort(j)} = K;
        end
    end   
end

if ~flag
    rez.R_CCG  = min(rez.R_CCG , rez.R_CCG');
    rez.Q_CCG  = min(rez.Q_CCG , rez.Q_CCG');
end