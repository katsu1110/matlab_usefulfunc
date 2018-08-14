function tu = encoding_tuning(stm, res)
% compute encoding properties that the tuning curve has
% assuming "rate encoding"
%
% INPUT:
% stm ... stimulus parameters
% res ... neural responses (e.g. firing rate)
% 
% OUTPUT:
% - tuning curve (average response to each stimulus)
% - reliability (var average / var all)
% - selectivity (max average res - min / mean)
% - discriminability (mutual information derived from decoding)
% - metabolic cost of encoding (response unpredictability)
%

ntr = size(stm, 1);
if ~isequal(ntr, size(res, 1))
    error('The size of stm and res must match!')
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
tu.selectivity = (max(tu.mean) - min(tu.mean))/mean([max(tu.mean), min(tu.mean)]);

%% 
% discriminability
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
tu.discriminability = 0;
for i = 1:lenuni
    for j = 1:lenuni
        Pij = sum(stm==i & pred==j)/ntr;
        Pi = sum(stm==i)/ntr;
        Pj = sum(pred==j)/ntr;
        if Pij > 0 && Pi > 0 && Pj > 0
            tu.discriminability = tu.discriminability ...
                + Pij*log2(Pij/(Pi*Pj));
        end
    end
end

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