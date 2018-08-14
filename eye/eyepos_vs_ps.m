function data = eyepos_vs_ps(x, y, p)
% check whether eye positions are affected by pupil size
% after Choe et al., 2016
% INPUT:
% x ... x position of the eye
% y ... y position of the eye
% p ... pupil size
%
% OUTPUT:
% data ... output structure
%

%%
% correlation
data.original = compute_correlation(x, y, p);

%%
% 2nd polynomial regression
data.polyreg_xparas = polyfit(p, x, 2);
data.correct_x = x - polyval(data.polyreg_xparas, p);
data.polyreg_yparas = polyfit(p, y, 2);
data.correct_y = y - polyval(data.polyreg_yparas, p);

%% 
% check no correlation between corrected eye pos and ps
data.corrected = compute_correlation(data.correct_x, data.correct_y, p);

function out = compute_correlation(x, y, p)
[r, pval] = corrcoef(x, p);
out.Pearson.r(1) = r(1,2);
out.Pearson.p(1) = pval(1,2);
[r, pval] = corrcoef(y, p);
out.Pearson.r(2) = r(1,2);
out.Pearson.p(2) = pval(1,2);
[out.Spearman.r(1), out.Spearman.p(1)]...
    = corr(x', p', 'type', 'Spearman');
[out.Spearman.r(2), out.Spearman.p(2)]...
    = corr(y', p', 'type', 'Spearman');