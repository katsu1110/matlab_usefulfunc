function unity_scatter(x, y, label)
%%
% generate a scatter plot with the unity line
% INPUT:
% x, y ... vector with the same length
% label ... interger vector with the same length, indicating a group 
%
% EXAMPLE: unity_scatter(randn(30,1), randn(30,1), randi(2, 30, 1)-1)
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% check legnth of x, y
if size(x, 1)==1; x = x'; end
if size(y, 1)==1; y = y'; end
if nargin < 2 && size(x, 2)==2; y = x(:, 2); x = x(:, 1); end
if nargin < 3; label = ones(size(x)); end
if length(x)~=length(y)
    error('vector x and y must be the same-length')
end
if length(label)~=length(x)
    error('label length must be the same as the data')
end

% range
all = [x(:); y(:)];
dist = max(all) - min(all);
range = [min(all) - 0.1*dist, max(all) + 0.1*dist];
% range(1) = floor(10*range(1))/10;
% range(2) = ceil(10*range(2))/10;

% unity
plot(range, [0, 0], '-','color',0.4*ones(1,3))
hold on;
plot([0, 0], range, '-','color',0.4*ones(1,3))
hold on;
plot(range, range, '-','color',0.4*ones(1,3))

% scatter
groups = unique(label);
n_group = length(groups);
map = [lines(n_group); [0 0 0]];
if n_group==1
    stats = cell(1, 1);
    map = map(end, :);
else
    stats = cell(1, n_group + 1);
end
for n = 1:n_group
    hold on;
    xd = x(label==groups(n)); yd = y(label==groups(n));
    [xd, yd] = nan_remove_pair(xd, yd);
    scatter(xd, yd, 20, 'filled','marker','o', 'markerfacecolor', map(n,:), ...
        'markerfacealpha',0.4, 'markeredgecolor','w','markeredgealpha',0.8)
    stats{n} = pair_tests([xd, yd]);    
end

% axis
axis([range range])
set(gca, 'box', 'off', 'TickDir', 'out')
set(gca, 'XTick', [range(1) range(end)])
set(gca, 'YTick', [range(1) range(end)])

% stats across groups
if n_group > 1
   hold on;
   xd = x; yd = y;
   [xd, yd] = nan_remove_pair(xd, yd);
   stats{n+1} = pair_tests([xd, yd]);
end

% display stats (non-parametric)
lens = length(stats);
for n = 1:lens
    if n==lens
        g_lab = 'all';
    else
        g_lab = ['group ' num2str(n)];
    end
    if stats{n}.pair(1).signrank.p < 0.05
        text(range(1)+0.025*(range(2)-range(1)), range(1)+(0.975 - 0.075*(n-1))*(range(2)-range(1)), ...
            [g_lab ': n=' num2str(stats{n}.pair(1).n) ', p<' num2str(pval_inequality(stats{n}.pair(1).signrank.p))], ...
            'fontsize', 6, 'color', map(n,:))
    else
        text(range(1)+0.025*(range(2)-range(1)), range(1)+(0.975 - 0.075*(n-1))*(range(2)-range(1)), ...
            [g_lab ': n=' num2str(stats{n}.pair(1).n)...
            ', p=' num2str(pval_inequality(stats{n}.pair(1).signrank.p))], 'fontsize', 6, 'color', map(n,:))
    end
end
        
% axis square

