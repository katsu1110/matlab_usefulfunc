function spk = mem2fire(varargin)
%% 
% membrane potential-based model, similar to Seillier et al., 2017
%

%% 
% stimulus type
stmtype = 'or';

% stimulus sequence
ntr = 100;
nframe = 1; 
stm_resolution = 10; % ms

% filter
filter_me = 57; % ms
filter_sd = 6;

% threshold for firing
c = 15;
v_thre = 0.05;

% visualize
plot_flag = 0;

j = 1;              
while j <= length(varargin)
    switch varargin{j}
        case 'stmtype'
            stmtype = varargin{j+1};            
        case 'ntr'
            ntr = varargin{j+1};
        case 'nframe'
            nframe = varargin{j+1};
        case 'stm_resolution'
            stm_resolution = varargin{j+1};
        case 'filter_me' 
            filter_me = varargin{j+1};
        case 'filter_sd'
            filter_sd = varargin{j+1};
        case 'c'
            c = varargin{j+1};       
        case 'v_thre'
            v_thre = varargin{j+1};    
        case 'plot'
            plot_flag = varargin{j+1};
    end
    j = j + 2;
end

% stimulus type and tuning in membrane potential
switch stmtype
    case 'or'
        SD = 40;
        AMP = 20;
        or = 22.5:22.5:180;
        v_tu = normpdf(or, 90, SD);
        v_tu = AMP*v_tu/max(v_tu);
    case 'co'
        N = 2;
        AMP = 20;
        c50 = 15;
        OFFSET = 3;
        or = linspace(1.56, 100, 7);
        v_tu = AMP*((or.^N)./(c50^N + or.^N)) + OFFSET;
end

%%
% stimulus sequence
stms = datasample(or, ntr*nframe, 'Replace', true)';
stm = zeros(ntr, nframe*stm_resolution);
begin = [1, 1];
for i = 1:nframe
   stm(:, begin(1):begin(1)+stm_resolution-1)= ...
        repmat(stms(begin(2):begin(2)+ntr-1), 1, stm_resolution);
   begin(1) = begin(1) + stm_resolution - 1;
   begin(2) = begin(2) + ntr - 1;
end

% stimulus-driven membrane potential
unistm = unique(or);
lenuni = length(unistm);
vm = stm;
for i = 1:lenuni
    vm(stm==or(i)) = v_tu(i);
end

% 1D Gaussian temporal filter
g_x = fspecial('gaussian',[1 filter_me], filter_sd);

% fluctuating membrane potential
vmf = vm;
for i = 1:ntr
    vmf(i,:) = imfilter(vm(i,:), g_x);
end

% spike rates
fr = c*(vmf - v_thre);
fr(fr <= 0) = 0;

% firing
spk = arrayfun(@(x) poissrnd(x), fr);

%%
% visualization
if ~ismember(plot_flag, 0)
    close all;
    figure;

    % stimulus sequence
    subplot(2,4,1)
    imagesc(stm)
    colorbar('southoutside')
    title('stimulus sequence')

    % stimulus tuning
    subplot(2,4,2)
    plot(or, v_tu, '-ok')
    title('Vm tuning')

    % temporal filter
    subplot(2,4,3)
    plot(g_x)
    title('temporal filter')

    % fluctuating membrane potential
    subplot(2,4,4)
    imagesc(vmf)
    colorbar('southoutside')
    title('membrane potentital')

    % spike rates
    subplot(2,4,5)
    imagesc(fr)
    colorbar('southoutside')
    title('spike rates')

    % firing
    subplot(2,4,6)
    imagesc(spk)
    colorbar('southoutside')
    title('firing')

    % tuning curve of firing rate
    subplot(2,4,7)
    tu = zeros(lenuni, 2);
    for i = 1:lenuni
        tu(i,1) = mean(mean(spk(stm(:,1)==unistm(i), :), 2));
        v = spk(stm(:,1)==unistm(i), :);
        tu(i,2) = std(v(:));
    end
    errorbar(or, tu(:,1), tu(:,2), '-ok', 'capsize', 0)
    title('spike tuning')
end