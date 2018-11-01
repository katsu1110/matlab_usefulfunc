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
% - discriminability (decoding accuracy)
% - metabolic cost of encoding (response unpredictability)
%

ntr = size(stm, 1);
if ~isequal(ntr, size(res, 1))
    error('The size of stm and res must match!')
end
if nargin < 3; stmtype = 'other'; end

if ~strcmp(stmtype, 'or')
    stmtype = 'other';
end

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
acc = zeros(ntr, 1);
tr = 1:ntr;
categ = categorical(stm);
for i = 1:ntr
    % leave-one-out
    trs = tr(~ismember(tr, i));

    % multinomial logistic regression
    B = mnrfit(res(trs), categ(trs));
    
    % model prediction of stimulus type
    prb = mnrval(B, res(i));
    acc(i) = ismember(find(tu.unistm==stm(i)), find(prb==max(prb)));
end
tu.discriminability = sum(acc)/ntr;

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
% metabolic cost (entropy, conditional entropy, mutual information)
tu.metabcost = zeros(1, 3);
int_res = round(res);
unires = unique(int_res);
lenr = length(unires);
for r = 1:lenr
    % entropy
    pr = sum(int_res==unires(r))/ntr;
    tu.metabcost(1) = tu.metabcost(1) + pr*log2(1/pr);
    
    % conditional entropy
    for s = 1:lenuni
        ps = sum(stm==tu.unistm(s))/ntr;
        prs = sum(int_res==unires(r) & stm==tu.unistm(s))/sum(stm==tu.unistm(s));
        if prs > 0
            tu.metabcost(2) = tu.metabcost(2) + ps*prs*log2(1/prs);
        end
    end
end
% mutual information
tu.metabcost(3) = tu.metabcost(1) - tu.metabcost(2);
 
% pred = zeros(ntr, 1);
% for i = 1:ntr
%     % leave-one-out
%     trs = tr(~ismember(tr, i));
%     
%     % train a optimal linear decoder (regression)
%     beta = glmfit(stm(trs), res(trs), 'normal', 'link', 'identity', 'constant', 'on');
%     
%     % trial-by-trial prediction
%     pred(i) = glmval(beta, stm(i), 'identity', 'constant', 'on');
% end
% tu.metabcost = var(abs(res - pred))/var(res);

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
