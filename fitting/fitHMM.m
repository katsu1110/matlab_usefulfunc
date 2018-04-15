function hmm_estimate = fitHMM(seq, n_state)
% fit the Hidden Markov Model (HMM) to sequence of (presumably) spike counts
% in order to estimate the HMM parameters
% INPUT: % INPUT: seq ... vector or matrix of spike counts 
%        (row: channel, column: spike-count x time)
%        ntr ... the number of trials
%        n_state ... the number of component
% OUTPUT: hmm_estimate ... output structure
%         n: the number of states, transition: transition matrix
%         emission: emission matrix, fr: firing rate in each state
%         likelihood: model likelihood, likelystates: estimated states
%         duration: duration of each state, variance_explained: variance
%         explained, err_...: confidence intervals of the HMM parameters
%         (this is not implemented...)
%
% NOTE: when the input seq is spike-count, it is highly recommended that
%       the raw spike-count is square-root-transformed to stabilize the variance
%
% written by Katsuhisa (15.04.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% (n_channel x frames --> 1 x frames) normalization
nc = size(seq, 1);
beta = ones(nc, 1);
if nc > 1
    for n = 1:nc
        beta(n) = mean(seq(n, :), 2);
        seq(n, :) = seq(n, :)./beta(n);
    end
    seq = mean(seq, 1);
end

% to avoid 0
seq = seq + 1;

% fit HMM 10 times to select sets of parameters yielding the largest likelihood 
% to avoid local optima
li = zeros(1, 10);
ttr_temp = cell(1, 10); emt_temp = cell(1, 10);
% ttr_ph = []; emt_ph = [];
for r = 1:10    
    % initial guess
    [tr_guess, em_guess] = params_initializer(seq, n_state);
    
    % fit HMM
    [ttr_temp{r}, emt_temp{r}] = hmmtrain(seq, tr_guess, em_guess, ...
        'Algorithm', 'BaumWelch', 'Maxiterations', 10000);  
%     ttr_ph = [ttr_ph; ttr_temp{r}(:)'];
%     emt_ph = [emt_ph; emt_temp{r}(:)'];
    
    % posterior probability
    posterior = hmmdecode(seq, ttr_temp{r}, emt_temp{r});
    li(r) = mean(abs(posterior(1,:) - 0.5)) + 0.5;
end
[li, maxidx] = max(li);
ttr = ttr_temp{maxidx};
emt = emt_temp{maxidx};
% ttr_low = reshape(quantile(ttr_ph, 0.05), n_state, []);
% ttr_high = reshape(quantile(ttr_ph, 0.95), n_state, []);
% emt_low = reshape(quantile(emt_ph, 0.05), n_state, []);
% emt_high = reshape(quantile(emt_ph, 0.95), n_state, []);

% estimate the states
likelystates = hmmviterbi(seq, ttr, emt);

% fr & duration of each state, variance explained
seq = seq - 1;
fr = zeros(nc, n_state);
frvec = seq;
for n = 1:n_state
    fr(:, n) = beta*mean(seq(likelystates==n));
    frvec(likelystates==n) = mean(fr(:, n), 1);
    dr.state(n).duration = [];
end
l = 1;
for i = 2:length(likelystates)
    if likelystates(i)==likelystates(i-1)
        l = l + 1;
    else
        dr.state(likelystates(i-1)).duration = ...
            [dr.state(likelystates(i-1)).duration, l];
        l = 1;
    end
end
 
% variance explained
var_tot = sum((seq - mean(seq)).^2);
var_res = abs(var_tot - sum((seq - frvec).^2));
varexp = 1 - (var_res/var_tot);
 
% % bootstrap x10 to estimate confidence intervals of the HMM parameters
% ttr_err = cell(1,10); emt_err = cell(1,10); fr_err = cell(1,10);
% for i = 1:10
%     start = randi(floor(ntr/2));
%     idx = start:start+floor(ntr/2)-1;
%     seq_temp = seq(idx);
% %     [ttr_err{i}, emt_err{i}] = hmmtrain(seq_temp, ttr, emt, ...
% %         'Algorithm', 'BaumWelch', 'Maxiterations', 1000);
%     [ttr_err{i}, emt_err{i}] = hmmestimate(seq_temp, likelystates(idx));
%     likelystates_temp = hmmviterbi(seq_temp, ttr_err{i}, emt_err{i});
%     seq_temp = seq_temp - 1;
%     fr_err_temp = ones(nc, n_state);
%     for n = 1:n_state
%         fr_err_temp(:, n) = beta*mean(seq_temp(likelystates_temp==n));
%     end
%     fr_err{i} = fr_err_temp;
% end
% % estimate confidence intervals
% mat_ttr = []; mat_emt = []; mat_fr = [];
% for i = 1:length(ttr_err)
%     mat_ttr = [mat_ttr, reshape(ttr_err{i}, [], 1)];
%     mat_emt = [mat_emt, reshape(emt_err{i}, [], 1)];
%     mat_fr = [mat_fr, reshape(fr_err{i}, [], 1)];
% end
% ttr_err{2} = reshape(quantile(mat_ttr', 0.95)', n_state, []);
% ttr_err{1} = reshape(quantile(mat_ttr', 0.05)', n_state, []);
% emt_err{2} = reshape(quantile(mat_emt', 0.95)', n_state, []);
% emt_err{1} = reshape(quantile(mat_emt', 0.05)', n_state, []);
% fr_err{2} = quantile(mat_fr, 0.95);
% fr_err{1} = quantile(mat_fr, 0.05);

% structurize
hmm_estimate.n_state = n_state;
hmm_estimate.trnsMat = ttr;
hmm_estimate.emtMat = emt;
hmm_estimate.fr = fr;
hmm_estimate.likelihood = li;
hmm_estimate.processed_seq = seq;
hmm_estimate.likelystates = likelystates;
hmm_estimate.duration = dr;
hmm_estimate.variance_explained = varexp;
% hmm_estimate.err(1).trnsMat = ttr_low;
% hmm_estimate.err(2).trnsMat = ttr_high;
% hmm_estimate.err(1).emtMat = emt_low;
% hmm_estimate.err(2).emtMat = emt_high;
% hmm_estimate.err.fr = fr_err;

% % subfunctions
% function [train, test] = seq_split(ntr, nframe, seed)
% rng(seed);
% train_tr = datasample(1:ntr,round(ntr/2),'Replace',false);
% train = []; test = [];
% for i = 1:ntr
%     if ismember(i, train_tr)
%         train = [train, i:(i+nframe-1)];
%     else
%         test = [test, i:(i+nframe-1)];
%     end
% end

function r = drchrnd(a,n)
% random sampling from dirichlet distribution
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

function [tr_guess, em_guess] = params_initializer(seq, n_state)
% uni = unique(seq);
uni = 1:max(seq);
tr_guess = []; em_guess = [];
for n = 1:n_state
    tr_guess = [tr_guess; drchrnd(100*ones(1,n_state)/n_state, 1)];
    em_guess = [em_guess; poisspdf(uni, quantile(seq, n_state/100))];
end