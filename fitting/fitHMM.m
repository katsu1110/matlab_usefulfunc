function [hmm_estimate] = fitHMM(spikecount, n_comp)
% fit the Hidden Markov Model (HMM) to sequence of spike counts
% in order to estimate the HMM parameters
% INPUT: spikecount ... vector or matrix of spike counts (integer)
%        n_comp ... the number of component (1, 2 or 3)
% OUTPUT: hmm_estimate ... output structure
%         n: the number of states, transition: transition matrix
%         emission: emission matrix, fr: firing rate in each state
%         likelihood: model likelihood, likelystates: estimated states
%         duration: duration of each state, variance_explained: variance
%         explained, err_...: confidence intervals of the HMM parameters
%
% written by Katsuhisa (23.02.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% initialization
likeli = 0;
ntr = size(spikecount,1);
idx = 1:ntr;

% to avoid 0
spikecount = spikecount + 1;

% one long vector
seq = mat2seq(spikecount);

% train, test 
rng(19891220);
train_idx = datasample(idx,round(ntr/2),'Replace',false);
test_idx = idx(~ismember(idx,train_idx));
seq_train = mat2seq(spikecount(train_idx,:));
seq_test = mat2seq(spikecount(test_idx,:));

% fit HMM 10 times to select sets of parameters yielding the largest likelihood 
% to avoid local optima
for r = 1:10    
    % initial guess
    [tr_guess, em_guess] = params_initializer(seq, n_comp);
    if r==1
        ttr = tr_guess;
        emt = em_guess;
    end
    
    % fit the data with 2-fold cross-validation  
    [ttr_temp1, emt_temp1] = hmmtrain(seq_train, ...
        tr_guess, em_guess, 'Algorithm', 'BaumWelch', 'Maxiterations', 10000);  
    [ttr_temp2, emt_temp2] = hmmtrain(seq_test, ...
        tr_guess, em_guess, 'Algorithm', 'BaumWelch', 'Maxiterations', 10000);
    
    % posterior probability
    posterior1 = hmmdecode(seq_test, ttr_temp1, emt_temp1);
    posterior2 = hmmdecode(seq_train, ttr_temp2, emt_temp2);
    posterior1 = mean(abs(posterior1(1,:)-0.5));
    posterior2 = mean(abs(posterior2(1,:)-0.5));
    
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
[dr] = state_dur(likelystates, n_comp);

% firing rate in each state
seq = seq - 1;
fr = zeros(size(emt,1),1);
for n = 1:n_comp
    fr(n) = mean(seq(likelystates==n));
end

% variance explained
varexp = varexp_comp(seq, likelystates, fr, n_comp);

% bootstrap x10 to estimate confidence intervals of the HMM parameters
ttr_err = cell(1,10); emt_err = cell(1,10); fr_err = cell(1,10);
for i = 1:10
    tr = datasample(idx,ntr,'Replace',true);
    seq_temp = mat2seq(spikecount(tr,:));
    [ttr_err{i}, emt_err{i}] = hmmtrain(seq_temp, ttr, emt, ...
        'Algorithm', 'BaumWelch', 'Maxiterations', 500);
    likelystates_temp = hmmviterbi(seq_temp, ttr_err{i}, emt_err{i});
    seq_temp = seq_temp - 1;
    fr_err_temp = zeros(size(emt_err{i},1),1);
    for n = 1:n_comp
        fr_err_temp(n) = mean(seq_temp(likelystates_temp==n));
    end
    fr_err{i} = fr_err_temp;
end
ttr_err = CIestimate(ttr_err);
emt_err = CIestimate(emt_err);
fr_err = CIestimate(fr_err);

% structurize
hmm_estimate = struct('n', n_comp, 'transition', ttr, 'emission', emt, 'fr', fr, ...
    'likelihood', likeli + 0.5, 'likelystates', likelystates, 'duration', dr, ...
    'variance_explained', varexp, 'err_transition', ttr_err, 'err_emission', emt_err, 'err_fr', fr_err);
hmm_estimate = hmm_estimate(1);

% subfunctions
function r = drchrnd(a,n)
% random sampling from dirichlet distribution
p = length(a);
r = gamrnd(repmat(a,n,1),1,n,p);
r = r ./ repmat(sum(r,2),1,p);

function [tr_guess, em_guess] = params_initializer(seq, n_comp)
% uni = unique(seq);
uni = 1:max(seq);
switch n_comp
    case 1
        tr_guess = 1;
        em_guess = poisspdf(uni, normrnd(mean(seq),std(seq)));
    case 2
        tr_guess = [drchrnd([85, 15], 1); drchrnd([15, 85], 1)];
        b = normrnd(0.3,0.1);
        while b <= 0
            b = normrnd(0.3,0.1);
        end
        em_guess = [poisspdf(uni, quantile(seq,b));
            poisspdf(uni, quantile(seq,1-b))];
    case 3
        tr_guess = [drchrnd([82, 9, 9], 1); drchrnd([9, 82, 9], 1); drchrnd([9, 9, 82], 1)]; 
        c = normrnd(0.25,0.1);
        d = normrnd(0.5,0.1);
        while c <= 0 || d <= 0
            c = normrnd(0.25,0.1);
            d = normrnd(0.5,0.1);
        end
        em_guess = [poisspdf(uni, quantile(seq,c));
            poisspdf(uni, quantile(seq,d));
            poisspdf(uni, quantile(seq,1-c-d))];
end

function [dr] = state_dur(states, n_comp)
% measure the duration of each state
for n = 1:n_comp
    dr.state(n).duration = [];
end
switch n_comp
    case 1
        dr.state(1).duration = length(states);
    case 2
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
        dr.state(1).duration = dr1; dr.state(2).duration = dr2;
    case 3
        if states(1) == 1
            dr1 = 1; dr2 = []; dr3 = [];
        elseif states(1) == 2
            dr1 = []; dr2 = 1; dr3 = [];
        elseif states(1) == 3
            dr1 = []; dr2 = []; dr3 = 1;
        end
        for i = 2:length(states)
            if states(i)==1 && states(i-1)==1
                dr1(end) = dr1(end) + 1;
            elseif states(i)==1 && states(i-1)~=1
                dr1 = [dr1, 1];
            elseif states(i)==2 && states(i-1)~=2
                dr2 = [dr2, 1];
            elseif states(i)==2 && states(i-1)==2
                dr2(end) = dr2(end) + 1;
            elseif states(i)==3 && states(i-1)==3
                dr3(end) = dr3(end) + 1;
            elseif states(i)==3 && states(i-1)~=3
                dr3 = [dr3, 1];
            end
        end
        dr.state(1).duration = dr1; dr.state(2).duration = dr2;
        dr.state(3).duration = dr3;
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
var_tot = sum((seq - mean(seq)).^2);
var_res = abs(var_tot - sum((seq - frvec).^2));
varexp = 1 - (var_res/var_tot);