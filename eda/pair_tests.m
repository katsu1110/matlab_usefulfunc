function stats = pair_tests(pair1, pair2)
%%
% compute statistics for pairwise comparision
% INPUT: pair1 ... N x 2 matrix (N is the number of observations)
%        if pair2 is also provided, 2-sample tests to see the mean and the
%        median difference between pair1 and 2 are performed
%
% OUTPUT: stats ... structure containing statistics
%
% EXAMPLE: stats = pair_tests(rand(30, 2), rand(40, 2));
%

if nargin < 2; pair2 = []; end

% remove nan
[X1, Y1] = nan_remove_pair(pair1(:,1), pair1(:,2));
stats.pair(1).n = length(X1);

% parametric test
[~, stats.pair(1).ttest.p, ~, stats.pair(1).ttest.stats]...
    = ttest(X1, Y1);

% non-parametric test
[stats.pair(1).signrank.p, ~, stats.pair(1).signrank.stats]...
    = signrank(X1, Y1);

% for pair 2, if exists
if ~isempty(pair2)
    % remove nan
    [X2, Y2] = nan_remove_pair(pair2(:,1), pair2(:,2));
    stats.pair(2).n = length(X2);

    % parametric test
    [~, stats.pair(2).ttest.p, ~, stats.pair(2).ttest.stats]...
        = ttest(X2, Y2);

    % non-parametric test
    [stats.pair(2).signrank.p, ~, stats.pair(2).signrank.stats]...
        = signrank(X2, Y2);
    
    % parametric 2-sample test
    [~, stats.ttest2.p, ~, stats.ttest2.stats]...
        = ttest2(Y1-X1, Y2-X2);
    
    % non-parametric 2-sample test
    [stats.ranksum.p, ~, stats.ranksum.stats]...
        = ranksum(Y1-X1, Y2-X2);
end