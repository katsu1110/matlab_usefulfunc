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
% - metabolic cost of encoding (entropy of the response)
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
for i = 1:ntr
    % leave-one-out
    trs = tr(~ismember(tr, i));
    % train a optimal linear decoder (one vs all)
    SvmModels = cell(1, lenuni);
    for u = 1:lenuni
        y = zeros(ntr-1, 1);
        y(stm(trs)==tu.unistm(u)) = 1;
        SvmModels{u} = fitcsvm(res(trs), y);
    end
    % SVM score in each model
    scores = zeros(lenuni, 1);
    for u = 1:lenuni
        [~, s] = predict(SvmModels{u}, res(i));
        scores(u) = abs(s(:,1));
    end    
    % model prediction of stimulus type
    [~, maxidx] = max(scores);
    pred(i) = tu.unistm(maxidx);
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
% metabolic cost
[~, binedges] = histcounts(res, lenuni);
binres = res;
for i = 1:ntr
    for u = 1:lenuni
        if res(i) >= binedges(u) && res(i) < binedges(u+1)
            binres(i) = u;
        end
    end
end
tu.metabcost = 0;
for u = 1:lenuni
    pu = sum(binres==u)/ntr;
    tu.metabcost = tu.metabcost - ...
        pu*log2(pu);
end
