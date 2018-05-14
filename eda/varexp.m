function variance_explained = varexp(yorig, ypred)
% compute variance explained
me = mean(yorig);
variance_explained = sum(sqrt((ypred - me).^2))/sum(sqrt((yorig - me).^2));