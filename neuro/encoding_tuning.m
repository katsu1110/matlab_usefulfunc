function tu = encoding_tuning(stm, res, stmtype)
% compute encoding properties that the tuning curve has
% assuming "rate encoding"
%
% INPUT:
% stm ... stimulus parameters
% res ... neural responses (e.g. firing rate)
% stmtype ... 'or', or 'other'
% 
% OUTPUT:
% - tuning curve (average response to each stimulus)
% - reliability (var average / var all)
% - selectivity (anova)
% - fano factor (var / mean)
% - discriminability (mutual information derived from decoding)
% - metabolic cost of encoding (response unpredictability)
%

ntr = size(stm, 1);
if ~isequal(ntr, size(res, 1))
    error('The size of stm and res must match!')
end
if nargin < 3; stmtype = 'other'; end

%%
% tuning curve
tu.unistm = unique(stm)';
lenuni = length(tu.unistm);
tu.mean = nan(1, lenuni);
tu.std = nan(1, lenuni);
tu.ntr = nan(1, lenuni);
for u = 1:lenuni
    tu.ntr(u) = sum(stm==tu.unistm(u));
    tu.mean(u) = mean(res(stm==tu.unistm(u)));
    tu.std(u) = std(res(stm==tu.unistm(u)));
end

%%
% reliability
tu.reliability = var(tu.mean)/var(res);

%%
% selectivity
tu.selectivity = anova1(res, stm,'off');

%%
% SNR2
tu.snr2 = (tu.mean./tu.std).^2;
tu.snr2(isnan(tu.snr2)) = 0;
tu.snr2(isinf(tu.snr2)) = 100;

%% 
% discriminability (an ideal observer's decoding)
pred = zeros(ntr, 1);
tr = 1:ntr;
t = templateSVM('Standardize', 1);
for i = 1:ntr
    % leave-one-out
    trs = tr(~ismember(tr, i));

    % multiclass svm classifier
    SvmModel = fitcecoc(res(trs), stm(trs), 'Learners', t, 'Coding', 'onevsall');
    
    % model prediction of stimulus type
    pred(i) = predict(SvmModel, res(i));
end
tu.discriminability = mean((pred - stm).^2)/std(stm);

% tu.discriminability = 0;
% for i = 1:lenuni
%     for j = 1:lenuni
%         Pij = sum(stm==i & pred==j)/ntr;
%         Pi = sum(stm==i)/ntr;
%         Pj = sum(pred==j)/ntr;
%         if Pij > 0 && Pi > 0 && Pj > 0
%             tu.discriminability = tu.discriminability ...
%                 + Pij*log2(Pij/(Pi*Pj));
%         end
%     end
% end

%%
% metabolic cost (response entropy = unpredictability)
pred = zeros(ntr, 1);
for i = 1:ntr
    % leave-one-out
    trs = tr(~ismember(tr, i));
    
    % train a optimal linear decoder (regression)
    beta = glmfit(stm(trs), res(trs), 'normal', 'link', 'identity', 'constant', 'on');
    
    % trial-by-trial prediction
    pred(i) = glmval(beta, stm(i), 'identity', 'constant', 'on');
end
tu.metabcost = var(abs(res - pred))/var(res);

%%
% stimulus specific quantity
switch stmtype
    case 'or'
        % circular variance --- Ringach et al. (2002)
        k = exp(1i*2*tu.unistm);
        R = sum(tu.mean.*k)/sum(tu.mean);
        tu.unique.circularvariance = 1 - abs(R);
        
        % direction selectivity
        [rp, irp] = max(tu.mean);
        pdeg = tu.unistm(irp); 
        if pdeg - 180 < 0
            udeg = pdeg + 180;
        else
            udeg = pdeg - 180;
        end
        ru = tu.mean(tu.unistm==udeg);
        tu.unique.directionsel = (rp - ru)/(rp + ru);
    otherwise
        tu.unique = nan;
end
