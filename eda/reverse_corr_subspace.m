function rcsub = reverse_corr_subspace(stmMat, actMat, wnd, stm_samprate, plot_flag)
% generic function to perform a reverse correlation subspace analysis
% INPUT: 
% stmMat ... row: trials, column: stimulus values (e.g. hdx, or...)
% actMat ... row: trials, column: numeric vector (e.g. spikes, LFP...)
% wnd ... analysis window (ms)
% stm_samprate ... sampling rate of stimulus presentation  
% plot_flag ... 0, no plot; 1, plot
%
% OUTPUT:
% rcsub ... activity elicited by stimulus
%
% written by Katsuhisa (06.04.18)
% ++++++++++++++++++++++++++++++++++++++++++++++++++
 
if nargin<2; error('Provide stmMat & actMat!'); end
if nargin<3; wnd = 150; end
if nargin<4; stm_samprate = 100; end
if nargin<5; plot_flag = 0; end

% experiment setting -------------------------------
% the number of trials
ntr = size(stmMat, 1);

% stimulus duration
nframes = size(stmMat, 2);
stmdur = 1000*nframes/stm_samprate;
ndata = size(actMat ,2);

% equalize sampling rate of stmMat & actMat
% (note that I assume the sampling rate for actMat >= stmMat)
% actMatBin = nan(ntr, nframes);
% nperbin = floor(size(actMat, 2)/nframes);
% begin = 1;
% for n = 1:nframes
%     actMatBin(:,n) = nanmean(actMat(:, begin:begin+nperbin-1), 2);
%     begin = begin + nperbin;
% end
upstmMat = zeros(ntr, ndata);
for n = 1:ntr
    upstmMat(n,:) = interp1(1:nframes, stmMat(n,:), ...
                linspace(1, nframes, ndata), 'nearest');
end

% reverse correlation --------------------
nframeperwnd = floor(ndata*wnd/stmdur);

% initialize output structure 
unique_stm = unique(round(1000*stmMat(:))/1000);
nstm = length(unique_stm);
rcsub = struct('stm', []);
for n = 1:nstm
    rcsub.stm(n).val = unique_stm(n);
    rcsub.stm(n).rcmat = [];
end
for i = 1:ntr
    stpos = 1;
    while stpos + nframeperwnd - 1 <= ndata
        stmidx = unique_stm == upstmMat(i, stpos);
        rcsub.stm(stmidx).rcmat = ...
            [rcsub.stm(stmidx).rcmat; actMat(i, stpos:stpos+nframeperwnd-1)];
        stpos = stpos + 1;
    end
end

% mean & SD of responses relative to the grand mean & SD
overallmean = zeros(nstm, nframeperwnd);
overallsd = zeros(nstm, nframeperwnd);
noverall = zeros(1, nstm);
for n = 1:nstm
    noverall(n) = size(rcsub.stm(n).rcmat, 1);
    overallmean(n, :) = noverall(n)*mean(rcsub.stm(n).rcmat, 1);
    overallsd(n, :) = noverall(n)*std(rcsub.stm(n).rcmat, [], 1);
end
overallmean = overallmean/sum(noverall);
overallsd = overallsd/sum(noverall);
grandmean = 1000*mean(overallmean(:));
grandsd = 1000*mean(overallsd(:));
for n = 1:nstm
    rcsub.stm(n).mean = 1000*overallmean(n, :);
    rcsub.stm(n).sd = 1000*overallsd(n,:);
     rcsub.stm(n).occurrence = noverall(n);
end

% visualization -----------------------------
if plot_flag==1
    close all;
    col = lines(nstm);
    h = figure;
    subplot(1,2,1)
    plot([0 nframeperwnd+1], grandmean.*[1 1], ':k')
    hold on;
    xlim([0 nframeperwnd+1])
    ylabel({'response', '(mean)'})
    xlabel('time (ms)')
    set(gca, 'XTick', [0 nframeperwnd+1], 'XTickLabel', [0 wnd])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    subplot(1,2,2)
    plot([0 nframeperwnd+1], grandsd.*[1 1], ':k')
    hold on;
    xlim([0 nframeperwnd+1])
    ylabel('(SD)')
    xlabel('time (ms)')
    set(gca, 'XTick', [0 nframeperwnd+1], 'XTickLabel', [0 wnd])
    set(gca, 'box', 'off'); set(gca, 'TickDir', 'out')
    pl = nan(1, nstm);
    leg = cell(1, nstm);
    for n = 1:nstm
        subplot(1,2,1)
        hold on;
        plot(1:nframeperwnd, rcsub.stm(n).mean, '-', 'color', col(n,:))
        subplot(1,2,2)
        hold on;
        pl(n) = plot(1:nframeperwnd, rcsub.stm(n).sd, '-', 'color', col(n,:));
        leg{n} = ['stm:' num2str(unique_stm(n))];
    end
    legend(pl, leg, 'location', 'northeast'); legend('boxoff')
end
    
