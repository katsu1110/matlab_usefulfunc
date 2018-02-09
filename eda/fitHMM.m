function [likelystates, ttr, emt] = fitHMM(spikecount, n_comp)
% fit the Hidden Markov Model (HMM) to sequence of spike counts
% in order to estimate the HMM parameters
% INPUT: spikecount ... vector or matrix of spike counts (integer)
%        n_comp ... the number of component (1, 2 or 3)
% OUTPUT: likelystates ... estimated states
%         ttr ... transition probability matrix
%         emt ... emission matrix
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% initial guess
switch n_comp
    case 1
        tr_guess = 1;
        em_guess = poisspdf(1:10, mean(spikecount(:)));
    case 2
        tr_guess = [0.9 0.1; 0.1 0.9];
        em_guess = [poisspdf(1:10, quantile(spikecount(:),0.3));
            poisspdf(1:10, quantile(spikecount(:),0.7))];
    case 3
        tr_guess = [0.8 0.1 0.1; 0.1 0.8 0.1; 0.1 0.1 0.8];
        em_guess = [poisspdf(1:10, quantile(spikecount(:),0.25));
            poisspdf(1:10, quantile(spikecount(:),0.5));
            poisspdf(1:10, quantile(spikecount(:),0.75))];
end

% fit the data
[ttr, emt] = hmmtrain(spikecount, tr_guess, em_guess, ...
    'Algorithm', 'Viterbi', 'Maxiterations', 500);

% estimate the states
likelystates = hmmviterbi(spikecount, ttr, emt);