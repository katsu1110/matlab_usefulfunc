function hmm_estimate = fitHMM(seq, ntr, n_state)
% fit the Hidden Markov Model (HMM) to sequence of (presumably) spike counts
% in order to estimate the HMM parameters
% INPUT: seq ... vector or matrix of spike counts (row: channel, column: integer vector)
%        ntr ... the number of trials
%        n_state ... the number of component
% OUTPUT: hmm_estimate ... output structure
%         n: the number of states, transition: transition matrix
%         emission: emission matrix, fr: firing rate in each state
%         likelihood: model likelihood, likelystates: estimated states
%         duration: duration of each state, variance_explained: variance
%         explained, err_...: confidence intervals of the HMM parameters
%
% NOTE: when the input seq is spike-count, it is highly recommended that
%       the raw spike-count is square-root-transformed to stabilize the variance
%
% written by Katsuhisa (14.04.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% initialization
likeli = 0;
nframe = floor(size(seq,2)/ntr);

% normalizing seq
nc = size(seq,1);
beta = nan(nc, 1);
for n = 1:nc
    beta(n) = mean(seq(n,:),2);
    seq(n, :) = seq(n, :)/beta(n);
end
seq = round(seq*min(beta));
beta = beta/min(beta);

% to avoid 0
seq = seq + 1;

% train, test 
[train, test] = seq_split(ntr, nframe, 19891220);

% fit HMM 10 times to select sets of parameters yielding the largest likelihood 
% to avoid local optima
for r = 1:10    
    % initial guess
    [tr_guess, em_guess] = params_initializer(seq, n_state);
    if r==1
        ttr = tr_guess;
        emt = em_guess;
    end
    
    try
        % fit the data with 2-fold cross-validation  
        [ttr_temp1, emt_temp1] = hmmtrain(seq(:, train), ...
            tr_guess, em_guess, 'Algorithm', 'BaumWelch', 'Maxiterations', 10000);  
        [ttr_temp2, emt_temp2] = hmmtrain(seq(:, test), ...
            tr_guess, em_guess, 'Algorithm', 'BaumWelch', 'Maxiterations', 10000);
    catch
        disp('error for fitting for the initialized parameters')
        continue
    end
    try
        % posterior probability
        posterior1 = hmmdecode(seq(:, test), ttr_temp1, emt_temp1);
        posterior2 = hmmdecode(seq(:, train), ttr_temp2, emt_temp2);
        posterior1 = mean(abs(posterior1(1,:)-0.5)) + 0.5;
        posterior2 = mean(abs(posterior2(1,:)-0.5)) + 0.5;
    catch
        disp('error for decoding for the estimated parameters ')
        continue
    end
    
    % weighted average
    ttr_temp = ttr_temp1*(posterior1/(posterior1+posterior2))...
        + ttr_temp2*(posterior2/(posterior1+posterior2));
    emt_temp = emt_temp1*(posterior1/(posterior1+posterior2))...
        + emt_temp2*(posterior2/(posterior1+posterior2));
    li = (posterior1 + posterior2)/2;
    
    % select sets of parameters when the likelihood is higher    
    if li > likeli
        likeli = li;
        ttr = ttr_temp;
        emt = emt_temp;
    end
end

% estimate the states
likelystates = hmmviterbi(seq, ttr, emt);

% fr & duration of each state, variance explained
seq = seq - 1;
fr = zeros(nc,n_state);
frvec = seq;
for n = 1:n_state
    fr(:, n) = beta.*mean(seq(:, likelystates==n), 2);
    frvec(:, likelystates==n) = fr(:, n);
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
var_tot = sum((seq(:) - mean(seq(:))).^2);
var_res = abs(var_tot - sum((seq(:) - frvec(:)).^2));
varexp = 1 - (var_res/var_tot);

% bootstrap x10 to estimate confidence intervals of the HMM parameters
ttr_err = cell(1,10); emt_err = cell(1,10); fr_err = cell(1,10);
for i = 1:10
    train = seq_split(ntr, nframe, i);
    seq_temp = seq(:, train);
    [ttr_err{i}, emt_err{i}] = hmmtrain(seq_temp, ttr, emt, ...
        'Algorithm', 'BaumWelch', 'Maxiterations', 1000);
    likelystates_temp = hmmviterbi(seq_temp, ttr_err{i}, emt_err{i});
    seq_temp = seq_temp - 1;
    fr_err_temp = seq_temp;
    for n = 1:n_state
        fr_err_temp(:, n) = beta.*mean(seq_temp(:, likelystates_temp==n), 2);
    end
    fr_err{i} = fr_err_temp;
end
% estimate confidence intervals
hb_ttr = ttr_err{1}; lb_ttr = hb_ttr;
hb_emt = emt_err{1}; lb_emt = hb_emt;
hb_fr = fr_err{1}; lb_fr = hb_fr;
mat_ttr = []; mat_emt = []; mat_fr = [];
for i = 1:length(ttr_err)
    mat_ttr = [mat_ttr, ttr_err{i}(:)];
    mat_emt = [mat_emt, emt_err{i}(:)];
    mat_fr = [mat_fr, fr_err{i}(:)];
end
for i = 1:length(ttr_err{1}(:))
    hb_ttr(i) = quantile(mat_ttr(i,:), 0.95);
    lb_ttr(i) = quantile(mat_ttr(i,:), 0.05);
    hb_emt(i) = quantile(mat_emt(i,:), 0.95);
    lb_emt(i) = quantile(mat_emt(i,:), 0.05);
    hb_fr(i) = quantile(mat_fr(i,:), 0.95);
    lb_fr(i) = quantile(mat_fr(i,:), 0.05);
end
ttr_err = {lb_ttr, hb_ttr};
emt_err = {lb_emt, hb_emt};
fr_err = {lb_fr, hb_fr};

% structurize
hmm_estimate = struct('n', n_state, 'transition', ttr, 'emission', emt, 'fr', fr, ...
    'likelihood', likeli, 'likelystates', likelystates, 'duration', dr, ...
    'variance_explained', varexp, 'err_transition', ttr_err, 'err_emission', emt_err, 'err_fr', fr_err);
hmm_estimate = hmm_estimate(1);

% subfunctions
function [train, test] = seq_split(ntr, nframe, seed)
rng(seed);
train_tr = datasample(1:ntr,round(ntr/2),'Replace',false);
train = []; test = [];
for i = 1:ntr
    if ismember(i, train_tr)
        train = [train, i:(i+nframe-1)];
    else
        test = [test, i:(i+nframe-1)];
    end
end

function r = drchrnd(a,n)
% random sampling from dirichlet distribution
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

function [tr_guess, em_guess] = params_initializer(seq, n_state)
% uni = unique(seq);
uni = 1:max(seq(:));
tr_guess = []; em_guess = [];
for n = 1:n_state
    tr_guess = [tr_guess; drchrnd(100*ones(1,n_state)/n_state, 1)];
    em_guess = [em_guess; poisspdf(uni, quantile(seq(1,:), n_state/100))];
end