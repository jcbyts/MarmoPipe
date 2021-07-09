function rez = find_merges_jly(rez, flag, debug)
% this function merges clusters based on template correlation

% this is a custom merge function. The basic principles are to merge units
% with similar waveforms, autocorrelograms, and refractory
% cross-correlograms
%
%
% This differs from the Kilosort2 default merge because it treats units
% that already have refractory-period violations separately from those that
% don't.
%
%


if nargin < 3
    debug = false; % debugging plots figures for each pair
end

ops = rez.ops;
dt = 1/1000;
noisecluster = 660; % we will be moving some clusters to noise / some spikes
ampthresh = 1; % if refractory period violations are this many sds away from the cluster mean, we will remove them

Xsim = rez.simScore; % this is the pairwise similarity score
Nk = size(Xsim,1); % number of clusters
Xsim = Xsim - diag(diag(Xsim)); % remove the diagonal of ones

%--- Step 1: clean up refractory period violations

% if a spike that violates a refractory period is more than one standard
% deviation from the cluster mean, get rid of it

unitIds = unique(rez.st3(:,2));
NC = numel(unitIds);
for cc = 1:NC
    
    s1ix = rez.st3(:,2)==unitIds(cc);
    sptimes = rez.st3(s1ix,1);
    spamps = rez.st3(s1ix,3);
    mamp = mean(spamps);
    sdev = std(spamps);
    
    % find all refractory violations
    rfviols = find(diff(sptimes)<dt);
    rempots = rfviols(:)+[1 0];
    nviol = numel(rfviols);
    
    ampdev = abs((spamps(rempots) - mamp) / sdev);
    if nviol==1
        ampdev = ampdev(:)';
    end
    [~, id] = max(ampdev, [], 2);
    
    maxmask = false(nviol, 2);
    ampmask = ampdev > ampthresh;
    
    maxmask(sub2ind([nviol, 2], 1:nviol, id(:)')) = true;
    
    toremove = rempots(ampmask & maxmask);
    
    rez.st3(toremove,2) = noisecluster;
   
end

% sort by firing rate first
nspk = zeros(Nk, 1);
for j = 1:Nk
    nspk(j) = sum(rez.st3(:,2)==j); % determine total number of spikes in each neuron
end
[~, isort] = sort(nspk); % we traverse the set of neurons in ascending order of firing rates

fprintf('initialized spike counts\n')

if ~flag
  % if the flag is off, then no merges are performed
  % this function is then just used to compute cross- and auto- correlograms
   rez.R_CCG = Inf * ones(Nk);
   rez.Q_CCG = Inf * ones(Nk);
   rez.K_CCG = {};
end

xc = rez.xcoords;
yc = rez.ycoords/3;

% --- Step 2: loop over units and perform merges
for j = 1:Nk
    if debug
        fprintf('Checking unit %d\n', isort(j))
    end
    
    s1ix = find(rez.st3(:,2)==isort(j));
    s1 = rez.st3(s1ix, 1)/ops.fs; % find all spikes from this cluster
    if numel(s1)~=nspk(isort(j))
        fprintf('lost track of spike counts') %this is a check for myself to make sure new cluster are combined correctly into bigger clusters
    end
    
    if flag && numel(s1) < 500 % can't do anything with only 500 spikes
        fprintf('setting %d to noise \n', isort(j))
        rez.st3(rez.st3(:,2)==isort(j), 2)=noisecluster;
        continue
    end
    
    % sort all the pairs of this neuron, discarding any that have fewer spikes
    [ccsort, ix] = sort(Xsim(isort(j),:) .* (nspk'>numel(s1)), 'descend');
    ienu = find(ccsort<.5, 1) - 1; % find the first pair which has too low of a correlation
   
    % extract waveform
    wf1 = rez.Wraw(:,:,isort(j));
    
    mwf1 = max(max(abs(wf1)));
    % --- THROW OUT UNITS WITH TINY WAVEOFMRS
    if mwf1 < 2 % hard-coded, might need to be a parameter
        fprintf('setting %d to noise because low amplitde\n', isort(j))
        rez.st3(rez.st3(:,2)==isort(j), 2)=noisecluster;
        continue
    end
    
    % find center channel
    amp1 = sqrt(sum(wf1.^2,2));
    [~,ctrCh] = max(amp1);
    
    chIdx = ctrCh + [-1 0 1];
    chIdx(chIdx > size(wf1, 1)) = [];
    chIdx(chIdx < 1) = [];
    
    % for all pairs above 0.5 correlation
    for k = 1:ienu
        s2ix = (rez.st3(:,2)==ix(k));
        if sum(s2ix)==0
            continue
        end
        s2 = rez.st3(s2ix, 1)/ops.fs; % find the spikes of the pair
        
        % extract waveform
        wf2 = rez.Wraw(:,:,ix(k));
        
        % center
        wc1 = wf1(chIdx,:);
        wc2 = wf2(chIdx,:);
        wfSSE = sum( (wc1(:)-wc2(:)).^2);
        wfR2 = 1 - wfSSE/sum( (wc1(:) - mean(wc1(:))).^2);
        wfCorr = corr(wf1(:), wf2(:));
        
        if debug
            figure(999); clf
            
            % plot waveforms across the probe
            subplot(1,2,1)
            
            % target neuron
            xax = linspace(-90, 90, size(rez.Wraw, 2));
            plot( (xax+xc)', wf1' - yc', 'b'); hold on
            xlabel('Time / Shank')
            ylabel('Amplitude / Depth')
            title('Templates')
            set(gca, 'XTickLabel', [], 'YTickLabel', [])
            axis tight
            
            % possible merge
            plot( (xax+xc)', wf2' - yc', 'r');
            
            subplot(3,2,2)
            plot(wf1(chIdx,:)' + (1:numel(chIdx))*50, 'b'); hold on
            plot(wf2(chIdx,:)' + (1:numel(chIdx))*50, 'r');
            title(sprintf('Centered: rho=%02.2f, r2=%02.2f', wfCorr, wfR2))
            set(gca, 'XTickLabel', [], 'YTickLabel', [])
            drawnow
        end
        
       
        mwf2 = max(max(abs(wf2)));
        fprintf('WF rho=%2.2f,\tr2=%2.2f\n', wfCorr, wfR2)
        fprintf('amp1=%2.2f,\tamp2=%2.2f\n', mwf1, mwf2)
        
        % if potential merge is low-amplitude, throw it out
        if mwf2 < 2
            fprintf('Setting low amplitude spike to noise\n')
            pause(0.1)
            rez.st3(s2ix, 2)=noisecluster;
            continue
        end
        
        if mwf2 < 3 && mwf1 < 3 % both are low-amplitude junk
            
            i = ix(k);
             
            % merge them
            collision = false(nspk(isort(j)),1); % who cares about detecting collisions for multi-units, they already both have major refractory period violations
            
            % now merge j into i and move on
            rez.st3(s1ix(~collision),2) = i; % simply overwrite all the spikes of neuron j with i (i>j by construction)
            rez.st3(s1ix(collision),2) = noisecluster; % set collisions to noise (will be removed)
            nspk(i) = nspk(i) + sum(~collision); % update number of spikes for cluster i
            
            fprintf('merged %d into %d \n', isort(j), i)
            
            break
        end
        
        if ((wfCorr > .9) || (wfR2 > .5))
            fprintf('Potential Merge...\n')
        else %don't even bother if the waveforms aren't remotely similar
            continue
        end
        
        
   
        % compute cross-correlograms, refractoriness scores (Qi and rir), and normalization for these scores
        nt = 50; % 50 ms lags
        Kac2 = ccg(s2, s2, nt, dt); % autocorrelation
        Kac1 = ccg(s1, s1, nt, dt);
        
        [K, Qi, Q00, Q01, rir] = ccg(s1, s2, nt, dt); % cross correlation
        
        % remove 0 lag from autocorrelation
        Kac1(nt+1) = nan;
        Kac2(nt+1) = nan;
        
        % normalize by edges of the autocorrelation (similar to poisosn
        % expectation)
        Kac1 = Kac1 / mean(Kac1([2 end-1]));
        Kac2 = Kac2 / mean(Kac2([2 end-1]));
        
        K = K / mean(K([2 end-1])); % normalize by edges
        
        
        % refractory period rate for ACGs and CCGs
        rf1 = Kac1(nt);
        rf2 = Kac2(nt);
        rfccg = K(nt);
        
        % is the ccg refractory
        if debug
            subplot(6,2,6)
            plot(-nt:nt, Kac1, 'b'); hold on
            plot(-nt:nt, Kac2, 'r');
            title('ACGs')
            
            subplot(6,2,8)
            plot(-nt:nt, K, 'k'); hold on
            if rfccg < .5
                plot(-1, K(nt), 'or')
            end
            title('CCG')
            
            subplot(3,2,6)
            plot( (Kac2 - Kac1).^2); hold on
            drawnow
            title([wfSSE nansum((Kac2 - Kac1).^2)])
        end
        
        % These are Kilosort2 metrics that they recommend
        Q = min(Qi/(max(Q00, Q01))); % normalize the central cross-correlogram bin by its shoulders OR by its mean firing rate
        R = min(rir); % R is the estimated probability that any of the center bins are refractory, and kicks in when there are very few spikes
        fprintf('Q=%2.2f, R=%2.2f\n', Q, R)
        
        fprintf('rf1=%2.2f, rf2=%2.2f, rfccg=%2.2f\n', rf1, rf2, rfccg)
        
%         input('Check')
        
        if flag
            
            if rf1 > .5 && rf2 > .5 % both multi unit (low threshold for merge)
                
                i = ix(k);
                
                if rfccg < 1 % merge any refractoriness in ccg
                    collision = false(nspk(isort(j)),1); % who cares about detecting collisions for multi-units, they already both have major refractory period violations
                    
                    % now merge j into i and move on
                    rez.st3(s1ix(~collision),2) = i; % simply overwrite all the spikes of neuron j with i (i>j by construction)
                    rez.st3(s1ix(collision),2) = noisecluster; % set collisions to noise (will be removed)
                    nspk(i) = nspk(i) + sum(~collision); % update number of spikes for cluster i
                    
                    fprintf('merged %d into %d \n', isort(j), i)
                    
                    break;
                end
                
            elseif (Q<1 && R<.05) || nansum( (K-Kac).^2 )<0.001 % if both refractory criteria are met
                
                i = ix(k);
                
                % find collisions: do not merge spikes that occured within
                % 1ms of eachother (they are the same spike)
                collision = false(nspk(isort(j)),1);
                for ispk = 1:nspk(isort(j))
                    delta = abs(s1(ispk) - s2) < 1e-3;
                    if  delta
                        collision(ispk) = true;
                    end
                end
                fprintf('%02.2f%% of spikes collided. Don''t merge those\n', 100*mean(collision))
                
                % now merge j into i and move on
                rez.st3(s1ix(~collision),2) = i; % simply overwrite all the spikes of neuron j with i (i>j by construction)
                rez.st3(s1ix(collision),2) = noisecluster; % set collisions to noise (will be removed)
                nspk(i) = nspk(i) + sum(~collision); % update number of spikes for cluster i
                
                fprintf('merged %d into %d \n', isort(j), i)
                
                break; % if a pair is found, we don't need to keep going (we'll revisit this cluster when we get to the merged cluster)
            end
        else
          % sometimes we just want to get the refractory scores and CCG
            rez.R_CCG(isort(j), ix(k)) = R;
            rez.Q_CCG(isort(j), ix(k)) = Q;

            rez.K_CCG{isort(j), ix(k)} = K;
            rez.K_CCG{ix(k), isort(j)} = K(end:-1:1); % the CCG is "antisymmetrical"
        end
    end
end

% clean up refractory period violations again

% if a spike that violates a refractory period is more than one standard
% deviation from the cluster mean, get rid of it

unitIds = unique(rez.st3(:,2));
NC = numel(unitIds);
for cc = 1:NC
    
    s1ix = rez.st3(:,2)==unitIds(cc);
    sptimes = rez.st3(s1ix,1);
    spamps = rez.st3(s1ix,3);
    mamp = mean(spamps);
    sdev = std(spamps);
    
    % find all refractory violations
    rfviols = find(diff(sptimes)<dt);
    rempots = rfviols(:)+[1 0];
    nviol = numel(rfviols);
    
    ampdev = abs((spamps(rempots) - mamp) / sdev);
    if nviol==1
        ampdev = ampdev(:)';
    end
    [~, id] = max(ampdev, [], 2);
    
    maxmask = false(nviol, 2);
    ampmask = ampdev > ampthresh;
    
    maxmask(sub2ind([nviol, 2], 1:nviol, id(:)')) = true;
    
    toremove = rempots(ampmask & maxmask);
    
    rez.st3(toremove,2) = noisecluster;
    
    
end

% remove noise spikes
noise = rez.st3(:,2)==noisecluster;
rez.cProj(noise,:) = [];
rez.cProjPC(noise,:,:) = [];
rez.st3(noise,:) = [];

if ~flag
    rez.R_CCG  = min(rez.R_CCG , rez.R_CCG'); % symmetrize the scores
    rez.Q_CCG  = min(rez.Q_CCG , rez.Q_CCG');
end