function fp = fixation_precision(x, y, p)
%% estimate fixation precision in the same way as Cherici et al.,2012
% ref: http://jov.arvojournals.org/article.aspx?articleid=2192348
% 
% INPUT:
% x: x-position of eye, y: y-position of eye, p: probability of 2D
% distribution of the eye
%
% OUTPUT: 
% fixspan ... output structure containing estimated fixation precision
%                  x: vector of eye-traces in x-pos
%                  y: vector of eye-traces in y-pos
%                  accuracy: fixation precision [deg^2]
%                  theta: degree to which a distribution of
%                  eye-position is tilted                        
%                  SI: symmetry index
%                  edge: for 2D histogram
% 
% NOTE: 
% - This function requires 'ndhist' 
% - x and y may be corrected before feeding into this function by substracting corresponding median
% - To visualize the 2D distribution, use 'sanePCcolor(x, y, p)'
%

% input
if nargin < 3; p = 0.75; end

% 2D histogram
[edgesX2,edgesY2,N] = ndhist(x, y);
close(gcf);
xrange = max(edgesX2) - min(edgesX2) + 1;
yrange = max(edgesY2) - min(edgesY2) + 1;

% convert counts into probability
Np = N/sum(N(:));

% cumulative probability
Np_des = sort(Np(:),'descend');
cum = zeros(1,length(Np_des)-1);
cum(1) = Np_des(1);
for t = 2:length(Np_des)-1
    cum(t) = Np_des(t) + cum(t-1);
end
idx = find(cum <= p, 1, 'first');
ratio = idx/(length(Np_des)-1);

% convert the computed ratio into view angle [degree^2]
p_area = ratio*xrange*yrange;
% disp(['fixation precision: ' num2str(p_area) ' degree^2'])

% PCA analysis ==============================
% covariance matrix
C = cov([x', y']);

% eigenvalue and eigenvector
[V, D] = eig(C);

% theta (degree to which the 2D distribution of the eye-positions tilts)
col = mod(find(D==max(D(:))), 2);
switch col
    case 0
            theta = atan2(V(2,2), V(1,2))*180/pi;
    case 1
            theta = atan2(V(2,1), V(1,1))*180/pi;
end
if theta < 0
        theta = theta + 180;
end

% disp(['theta: ' num2str(theta)])

% symmetry index (degree to which x and y eye-positions are symmetric)
ds = sort(D(:),'descend');
SI = sqrt(ds(2)/ds(1));

% disp(['symmetry index: ' num2str(SI)])

% structurize
fp = struct('accuracy', p_area, 'variance', [var(x), var(y)],'theta', theta, 'SI', SI, ...
    'edge_x', edgesX2, 'edge_y', edgesY2, 'Prob', Np);