function unity_scatter(x,y)
%%
% generate a scatter plot with the unity line
% INPUT:
% x, y ... vector or matrix (obs x feature)
% params ... cell array containing figure parameters
%
% example: unity_scatter(randn(30,3),randn(30,3))
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% check legnth of x, y
if size(x, 1) ~= size(y, 1) || size(x, 2) ~= size(y, 2) 
    error('vector x and y must be the same-length')
end

% range
all = [x(:); y(:)];
dist = max(all) - min(all);
range = [min(all) - 0.05*dist, max(all) + 0.05*dist];

% unity
plot(range, [0, 0], '-','color',0.4*ones(1,3))
hold on;
plot([0, 0], range, '-','color',0.4*ones(1,3))
hold on;
plot(range, range, '-','color',0.4*ones(1,3))

% scatter
dist = range(end) - range(1);
nobs = size(x, 2);
map = hsv(nobs);
for i = 1:nobs
    hold on;
    scatter(x(:,i), y(:,i), 50, 'filled','marker','o', 'markerfacecolor', map(i,:), ...
        'markerfacealpha',0.4, 'markeredgecolor','w','markeredgealpha',0.8)
    p = signrank(x(:,i), y(:,i));
    if p < 0.05
        text(0.65*dist+range(1), 0.05*i*dist+range(1),['p' num2str(i) '< ' num2str(pval_inequality(p))])
    else
        text(0.65*dist+range(1), 0.05*i*dist+range(1),['p' num2str(i) '=' num2str(pval_inequality(p))])
    end
end
text(0.65*dist+range(1),0.05*(i+1)*dist+range(1),['n=' num2str(length(x))])

% axis
axis([range range])
set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
set(gca, 'XTick', [range(1) range(end)])
set(gca, 'YTick', [range(1) range(end)])
% axis square
