function hmm_estimate = fitHMM(seq, n_state, analysis_channels)
% fit the Hidden Markov Model (HMM) to sequence of (presumably) spike counts
% in order to estimate the HMM parameters
%
% INPUT: seq ... matrix of spike counts (trials x spike counts x channels) 
%        n_state ... the number of component (2, in default)
%        analysis_channels ... channel number used for parameter estimation
% 
% OUTPUT: hmm_estimate ... output structure
%         n: the number of states, transition: transition matrix
%         emission: emission matrix, fr: firing rate in each state
%         likelihood: model likelihood, likelystates: estimated states
%         duration: duration of each state, variance_explained: variance
%         explained, err_...: confidence intervals of the HMM parameters
%         (this is not implemented...)
%
% % NOTE: - The algorithm works well if seq is binary (0 or 1).
%         - To transform the dimension of likelystates back, use
%         "reshape(likelystates, length(spike counts), length(trials))'"
%
% written by Katsuhisa (29.09.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% inputs
if nargin < 2; n_state = 2; end
if nargin < 3; analysis_channels = 1:size(seq, 3); end

% spike counts across channels
sc = zeros(size(seq, 1), size(seq, 2));
nc = length(analysis_channels);
for n = 1:nc
    sc = sc + seq(:,:,analysis_channels(n));
end

% 1d sequence
sc = sc';
sc = sc(:)';

% avoid 0, just for matlab
sc = sc + 1;

% fit HMM 10 times to select sets of parameters yielding the largest likelihood 
% to avoid local optima
repeat = 10;
li = zeros(1, repeat);
ttr_temp = cell(1, repeat); emt_temp = cell(1, repeat);
for r = 1:10    
    % initial guess
    [tr_guess, em_guess] = params_initializer(sc, n_state);
    
    % fit HMM
    [ttr_temp{r}, emt_temp{r}] = hmmtrain(sc, tr_guess, em_guess, ...
        'Algorithm', 'BaumWelch', 'Maxiterations', 10000);
    
    % posterior probability
    posterior = hmmdecode(sc, ttr_temp{r}, emt_temp{r});
    li(r) = mean(abs(posterior(1,:) - 0.5)) + 0.5;
end
[li, maxidx] = max(li);
ttr = ttr_temp{maxidx};
emt = emt_temp{maxidx};

% estimate the states
likelystates = hmmviterbi(sc, ttr, emt);

% fr & duration of each state, variance explained
frvec = sc - 1;
for n = 1:n_state
    hmm_estimate.state(n).fr = mean(sc(likelystates==n) - 1);
    frvec(likelystates==n) = hmm_estimate.state(n).fr;
    hmm_estimate.state(n).duration = [];
end
l = 1;
for i = 2:length(likelystates)
    if likelystates(i)==likelystates(i-1)
        l = l + 1;
    else
        hmm_estimate.state(likelystates(i-1)).duration = ...
            [hmm_estimate.state(likelystates(i-1)).duration, l];
        l = 1;
    end
end
 
% variance explained
sc = sc - 1;
var_tot = sum((sc - mean(sc)).^2);
var_res = abs(var_tot - sum((sc - frvec).^2));
hmm_estimate.variance_explained = 1 - (var_res/var_tot);
 
% structurize
hmm_estimate.n_state = n_state;
hmm_estimate.trnsMat = ttr;
hmm_estimate.emtMat = emt;
hmm_estimate.likelihood = li;
hmm_estimate.processed_seq = sc;
hmm_estimate.likelystates = likelystates;

function r = drchrnd(a,n)
% random sampling from dirichlet distribution
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

function [tr_guess, em_guess] = params_initializer(sc, n_state)
% initialize parameters
uni = 1:max(sc);
tr_guess = []; em_guess = [];
for n = 1:n_state
    tr_guess = [tr_guess; drchrnd(100*ones(1,n_state)/n_state, 1)];
    em_guess = [em_guess; poisspdf(uni, quantile(sc, n_state/100))];
end