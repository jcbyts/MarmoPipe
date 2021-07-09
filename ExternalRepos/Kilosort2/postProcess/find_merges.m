function rez = find_merges(rez, flag)

ops = rez.ops;
dt = 1/1000;

Xsim = rez.simScore;
Nk = size(Xsim,1);
Xsim = Xsim - diag(diag(Xsim));

% sort by firing rate first
nspk = zeros(Nk, 1);
for j = 1:Nk
    nspk(j) = sum(rez.st3(:,2)==j);        
end
[~, isort] = sort(nspk);

fprintf('initialized spike counts\n')

if ~flag
   rez.R_CCG = Inf * ones(Nk);
   rez.Q_CCG = Inf * ones(Nk);
   rez.K_CCG = {};
end

for j = 1:Nk
    s1 = rez.st3(rez.st3(:,2)==isort(j), 1)/ops.fs;
    if numel(s1)~=nspk(isort(j))
        fprintf('lost track of spike counts')
    end    
    [ccsort, ix] = sort(Xsim(isort(j),:) .* (nspk'>numel(s1)), 'descend');
    ienu = find(ccsort<.5, 1) - 1;
    
    for k = 1:ienu
        s2 = rez.st3(rez.st3(:,2)==ix(k), 1)/ops.fs;
        [K, Qi, Q00, Q01, rir] = ccg(s1, s2, 500, dt);
        Q = min(Qi/(max(Q00, Q01)));
%         Q = min(Qi/Q01);
        R = min(rir);
        
        if flag
            if Q<.2 && R<.05
                i = ix(k);
                % now merge j into i and move on
                rez.st3(rez.st3(:,2)==isort(j),2) = i;
                nspk(i) = nspk(i) + nspk(isort(j));
                fprintf('merged %d into %d \n', isort(j), i)
                % YOU REALLY SHOULD MAKE SURE THE PC CHANNELS MATCH HERE
                break;
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