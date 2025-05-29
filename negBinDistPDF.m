function p=negBinDistPDF(x,meanVal,varianceVal)
%NEGBINDISTPDF returns probability density values for a negative binomial distribution,
%parametrized by its mean and variance. This is a
%different parametrization than the standard one, but it is
%handy if the mean and variance are easily calculated.
% I took the conversion from mean and variance to r and p from:
% https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4

r=meanVal.*meanVal./(varianceVal-meanVal);  % transform to the standard r,p parametrization
% p=1 - (meanVal./varianceVal); % this was the equation in the genomebiology paper, but it is incorrect.
p=meanVal./varianceVal; % this was the equation from Wikipedia, which works

%p=(r+x-1 choose x) * p^r * (1-p)^x;  OR with the continuous gamma function:
%p=gamma(r+x)/gamma(r)/gamma(x+1) * p^r * (1-p)^x
logp= gammaln(r+x) - gammaln(r) - gammaln(x+1) + r.*log(p) + x.*log(1-p);
p=exp(logp);