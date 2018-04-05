function interaction_plot(mat)
% scatter plot to visualize an interaction effect
% INPUT: mat ... matrix where row represents observations and
%        column must be four-dimension vector (column 1 - 4)
%        The relationship of the columns must be summarized as following:
%               group A     group B
% effect 1        1           3
% effect 2        2           4
%
% This function creates a scatter where x-axis is effect 1 and y-axis is
% effect 2.
%
% written by Katsuhisa (05.04.18)
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

params0 = {30, 'filled','marker','o','markerfacecolor','k','markerfacealpha',0.8,...
    'markeredgecolor','w','markeredgealpha',0.8};
params1 = {30, 'filled','marker','^','markerfacecolor','k','markerfacealpha',0.8,...
    'markeredgecolor','w','markeredgealpha',0.8};

% range
minval = min(mat(:));
maxval = max(mat(:));
dist = maxval - minval;
range = [minval - 0.05*dist, maxval + 0.05*dist];

% unity
plot(range, [0, 0], '-','color',0.4*ones(1,3))
hold on;
plot([0, 0], range, '-','color',0.4*ones(1,3))
hold on;
plot(range, range, '-','color',0.4*ones(1,3))
hold on;

% scatter
nob = size(mat, 1);
for i = 1:nob
    plot(mat(i,[1 3]), mat(i,[2 4]), '-k')
    hold on;
    scatter(mat(i,1), mat(i,2), params0{:})
    hold on;
    scatter(mat(i,3), mat(i,4), params1{:})
    hold on;
end

% axis
axis([range range])
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
set(gca, 'XTick', [range(1) range(end)])
set(gca, 'YTick', [range(1) range(end)])
% axis square

% stats
dist = range(end) - range(1);
% [~,p] = ttest(x,y);
text(0.65*dist+range(1),0.1*dist+range(1),['n=' num2str(nob)])
p = signrank(mat(:,1),mat(:,2));
text(0.65*dist+range(1),0.2*dist+range(1),['p1=' num2str(pval_inequality(p))])
p = signrank(mat(:,3),mat(:,4));
text(0.65*dist+range(1),0.3*dist+range(1),['p2=' num2str(pval_inequality(p))])
