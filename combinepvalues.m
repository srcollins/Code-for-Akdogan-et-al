function p=combinepvalues(pvals,scale)
%uses a formula to combine a set of pvalues into one p-value. It can accept
%(and correspondingly return) regular or log p-values, as indicated by the
%scale variable. The formula used here was obtained from:
%http://www.loujost.com/Statistics%20and%20Physics/Significance%20Levels/CombiningPValues.htm
%the algorithm was obtained from:
%http://bioinformatics.oxfordjournals.org/cgi/reprint/14/1/48

if nargin<2
    scale='linear';
end

n=length(pvals(:));
    if strcmp(scale,'linear') & min(pvals)==0
        p=0;
    else
        if strcmp(scale,'log')
            logprod=sum(pvals(:));
        else
            logprod=sum(log(pvals(:)));
        end
        prod=exp(logprod);
        if n>1
            x=-1*logprod;
        end
        t=prod;
        q=prod;
        for i=1:(n-1)
            t=t*x/i;
            q=q+t;
        end
    end
    if strcmp(scale,'log')
        p=log(q);
    else
        p=q;
    end
