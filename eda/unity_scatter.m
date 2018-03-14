function unity_scatter(x,y,params)
%%
% generate a scatter plot with the unity line
% INPUT:
% x, y ... vectors with the same length
% params ... cell array containing figure parameters
%
% example: unity_scatter(randn(1,3),randn(1,3),{'markerfacecolor','r','markerfacealpha',0.4})
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% check legnth of x, y
if length(x) ~= length(y)
    error('vector x and y must be the same-length vectors')
end

% input
if size(x,1) > 1
    x = x';
end
if size(y,1) > 1
    y = y';
end
if nargin < 3
    params = {50, 'filled','marker','o','markerfacecolor','r','markerfacealpha',0.4,...
        'markeredgecolor','w','markeredgealpha',0.8};
end

% range
dist = max([x y]) - min([x y]);
range = [min([x y]) - 0.05*dist, max([x y]) + 0.05*dist];

% unity
plot(range, range, '-','color',0.4*ones(1,3))
hold on;

% scatter
scatter(x, y, params{:})

% axis
axis([range range])
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
set(gca, 'XTick', [range(1) range(end)])
set(gca, 'YTick', [range(1) range(end)])
% axis square

% stats
dist = range(end) - range(1);
% [~,p] = ttest(x,y);
text(0.65*dist+range(1),0.1*dist+range(1),['n=' num2str(length(x))])
p = signrank(x,y);
text(0.65*dist+range(1),0.2*dist+range(1),['p=' num2str(p)])
