function cm = confusion_matrix(vecs, corrtype)
% generate confusion matrix (correlation coefficient between pairs of the
% given vectors)
% INPUT:
% vecs ... cell arrray each containing vector (vectors have to be the same
% length)
% corrtype ... 'Pearson' or 'Spearman'
%

lenv = length(vecs);
cm = nan(lenv, lenv);
