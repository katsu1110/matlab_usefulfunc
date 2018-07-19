function cm = confusion_matrix(vecs, corrtype, mapcolor)
% generate confusion matrix (correlation coefficient between pairs of the
% given vectors)
% INPUT:
% vecs ... cell arrray each containing vector (vectors have to be the same
% length and vertical)
% corrtype ... 'Pearson' or 'Spearman'
% mapcolor ... color scheme of confusion matrix ('parula' in default)
%

if nargin < 2
    corrtype = 'Pearson';
end
if nargin < 3
    mapcolor = 'parula';
end
lenv = length(vecs);
cm = nan(lenv, lenv);
for r = 1:lenv
    for c = 1:lenv
        cm(r,c) = corr(vecs{r}, vecs{c}, 'type', corrtype);
    end
end
imagesc(1:lenv, 1:lenv, cm)
colormap(mapcolor)
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')