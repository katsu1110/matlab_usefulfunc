function [binned] = tcbin(v, binsize)
%%
% bin the vector v based on the specified binsize.
% note that this binning simply splits 'v' into 'binsize' elements
% from the beginning to the end (different from a 'hist' like method).
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++

lenv = length(v);
r = mod(lenv, binsize);
binned = floor(lenv/binsize)*ones(1, binsize);
i = 1;
while r > 0
    binned(i) = binned(i) + 1;
    r = r - 1;
    if i <= binsize
        i = i + 1;
    else
        i = 1;
    end
end
begin = 1;
for b = 1:binsize
    step = binned(b);
    binned(b) = nanmean(v(begin:begin+step-1));
    begin = begin + step;
end


