function [x, y] = nan_remove_pair(x, y, replace)
%%
% remove nans from paired-data
%
% INPUT: x, y ... vectors with the same length
%        replace ... 'none' (remove nans)
%                    'mean' (replace nans with mean)
%                    'median' (replace nans with median)
%

if nargin < 3; replace = 'none'; end

nans = isnan(x) | isnan(y);

switch replace
    case 'none'
        x(nans) = []; y(nans) = [];
    case 'mean'
        x(nans) = nanmean(x); 
        y(nans) = nanmean(y);
    case 'median'
        x(nans) = nanmedian(x); 
        y(nans) = nanmedian(y);
end