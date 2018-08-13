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
r = corrcoef(x, p);
data.Pearson(1) = r(1,2);
r = corrcoef(y, p);
data.Pearson(2) = r(1,2);

%%
% 2nd polynomial regression
data.polyreg_xparas = polyfit(p, x, 2);
data.correct_x = x - polyval(data.polyreg_xparas, p);
data.polyreg_yparas = polyfit(p, y, 2);
data.correct_y = y - polyval(data.polyreg_yparas, p);