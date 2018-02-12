function [hmm_estimate] = fitHMM(spikecount, n_comp)
% fit the Hidden Markov Model (HMM) to sequence of spike counts
% in order to estimate the HMM parameters
% INPUT: spikecount ... vector or matrix of spike counts (integer)
%        n_comp ... the number of component (1, 2 or 3)
% OUTPUT: likelystates ... estimated states
%         ttr ... transition probability matrix
%         fr ... firing rate in each state
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% fit HMM 10 times to select sets of parameters yielding the largest likelihood 
% to avoid local optima
likeli = 0;
ntr = size(spikecount,1);
idx = 1:ntr;
spikecount = spikecount + 1;
seq = mat2seq(spikecount);
uni = unique(seq);
for r = 1:10
    % random initial guess
    switch n_comp
        case 1
            tr_guess = 1;
            em_guess = poisspdf(uni, normrnd(mean(seq),std(seq)));
        case 2
            a = normrnd(0.9,0.05,2,1);
            tr_guess = [a(1), 1-a(1); 1-a(2), a(2)];
            b = normrnd(0.3,0.1);
            em_guess = [poisspdf(uni, quantile(seq,b));
                poisspdf(uni, quantile(seq,1-b))];
        case 3
            a = normrnd(0.8,0.05,3,1);
            b = normrnd(0.1,0.02,3,1);
            tr_guess = [a(1) b(1) 1-a(1)-b(1); 1-a(2)-b(2) a(2) b(2); b(3) 1-a(3)-b(3) a(3)];
            c = normrnd(0.25,0.1);
            d = normrnd(0.5,0.1);
            em_guess = [poisspdf(uni, quantile(seq,c));
                poisspdf(uni, quantile(seq,d));
                poisspdf(uni, quantile(seq,1-c-d))];
    end
    
    % fit the data with 2-fold cross-validation    
    train_idx = datasample(idx,round(ntr/2),'Replace',false);
    test_idx = idx(~ismember(idx,train_idx));
    [ttr_temp1, emt_temp1] = hmmtrain(mat2seq(spikecount(train_idx,:)), ...
        tr_guess, em_guess, 'Algorithm', 'BaumWelch', 'Maxiterations', 10000);    
    [ttr_temp2, emt_temp2] = hmmtrain(mat2seq(spikecount(test_idx,:)), ...
        tr_guess, em_guess, 'Algorithm', 'BaumWelch', 'Maxiterations', 10000);    
    ttr_temp = (ttr_temp1 + ttr_temp2)/2;
    emt_temp = (emt_temp1 + emt_temp2)/2;
    
    % posterior probability
    posterior = hmmdecode(seq, ttr_temp, emt_temp);
    li = mean(abs(posterior(1,:)-0.5));
    
    % select sets of parameters when the log-likelihood is high
    if li > likeli
        likeli = li;
        ttr = ttr_temp;
        emt = emt_temp;
    end
end

% estimate the states
likelystates = hmmviterbi(seq, ttr, emt);
[dr1, dr2] = state_dur(likelystates);

% firing rate in each state
fr = emt2fr(emt);

% variance explained
seq = seq - 1;
varexp = varexp_comp(seq, likelystates, fr, n_comp);

% bootstrap x10 to estimate confidence intervals of the HMM parameters
ttr_err = cell(1,10); emt_err = cell(1,10); fr_err = cell(1,10);
for i = 1:10
    tr = datasample(idx,ntr,'Replace',true);
    [ttr_err{i}, emt_err{i}] = hmmtrain(mat2seq(spikecount(tr,:)), ttr, emt, ...
        'Algorithm', 'Viterbi', 'Maxiterations', 500);
    fr_err{i} = emt2fr(emt_err{i});
end
ttr_err = CIestimate(ttr_err);
emt_err = CIestimate(emt_err);
fr_err = CIestimate(fr_err);

% structurize
hmm_estimate = struct('n', n_comp, 'transition', ttr, 'emission', emt, 'fr', fr, ...
    'loglikelihood', likeli + 0.5, 'likelystates', likelystates, 'duration0', dr1, 'duration1', dr2, ...
    'variance_explained', varexp, 'err_transition', ttr_err, 'err_emission', emt_err, 'err_fr', fr_err);
hmm_estimate = hmm_estimate(1);

% subfunctions
function fr = emt2fr(emt)
% convert emission matrix to firing rate (spike count/10ms)
fr = zeros(size(emt,1),1);
for i = 1:size(size(emt,2))
    fr(:,i) = fr(:,i) + (i - 1)*emt(:,i);
end

function [dr1, dr2] = state_dur(states)
% measure the duration of each state
if states(1) == 1
    dr1 = 1; dr2 = [];
elseif states(1) == 2
    dr1 = []; dr2 = 1;
end
for i = 2:length(states)
    if states(i)==1 && states(i-1)==1
        dr1(end) = dr1(end) + 1;
    elseif states(i)==1 && states(i-1)==2
        dr1 = [dr1, 1];
    elseif states(i)==2 && states(i-1)==1
        dr2 = [dr2, 1];
    elseif states(i)==2 && states(i-1)==2
        dr2(end) = dr2(end) + 1;
    end
end

function ci = CIestimate(cellmat)
% estimate confidence intervals
hb = cellmat{1}; lb = hb;
mat = [];
for i = 1:length(cellmat)
    mat = [mat, cellmat{i}(:)];
end
for i = 1:length(cellmat{1}(:))
    hb(i) = quantile(mat(i,:), 0.95);
    lb(i) = quantile(mat(i,:), 0.05);
end
ci = {lb, hb};

function seq = mat2seq(mat)
% tr x time --> one vector
mat = mat'; seq = mat(:)';

function varexp = varexp_comp(seq, likelystates, fr, n_comp)
% variance explained of the HMM
frvec = likelystates;
for n = 1:n_comp
    frvec(likelystates==n) = fr(n);
end
var_res = sum((seq - frvec).^2);
var_tot = sum((seq - mean(seq)).^2);
varexp = 1 - (var_res/var_tot);