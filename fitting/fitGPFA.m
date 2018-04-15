function gpr_estimate = fitGPFA(seq)
% fit the Gaussian Process Factor Analysis to sequence of (presumably) spike counts
% in order to estimate the HMM parameters
% INPUT: seq ... vector or matrix of spike counts 
%        (row: channel, column: spike-count x time)
% OUTPUT: gpr_estimate ... output structure
%
% NOTE: when the input seq is spike-count, it is highly recommended that
%       the raw spike-count is square-root-transformed to stabilize the variance 
% 
% written by Katsuhisa (15.04.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% (n_channel x frames --> 1 x frames) normalization
if size(seq, 1) > 1
    for n = 1:size(seq, 1)
        seq(n, :) = seq(n, :)./mean(seq(n, :), 2);
    end
    seq = mean(seq, 1);
end

% fit GPFA  
gprMdl = fitrgp([1:length(seq)]', seq', 'Basis', 'linear', 'KernelFunction', 'squaredexponential');

% loss of prediction
L = loss(gprMdl, [1:length(seq)]', seq');

% prediction
predseq = resubPredict(gprMdl);

% structurize
gpr_estimate.Mdl = gprMdl;
gpr_estimate.loss = L;
gpr_estimate.processed_seq = seq;
gpr_estimate.pred = predseq;