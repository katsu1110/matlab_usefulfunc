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
[r, pval] = corrcoef(x, p);
data.Pearson.r(1) = r(1,2);
data.Pearson.p(1) = pval(1,2);
[r, pval] = corrcoef(y, p);
data.Pearson.r(2) = r(1,2);
data.Pearson.p(2) = pval(1,2);
[data.Spearman.r(1), data.Spearman.p(1)]...
    = corr(x', p', 'type', 'Spearman');
[data.Spearman.r(2), data.Spearman.p(2)]...
    = corr(y', p', 'type', 'Spearman');

%%
% 2nd polynomial regression
data.polyreg_xparas = polyfit(p, x, 2);
data.correct_x = x - polyval(data.polyreg_xparas, p);
data.polyreg_yparas = polyfit(p, y, 2);
data.correct_y = y - polyval(data.polyreg_yparas, p);