function varModels= CLO_compute_variance_models(d,samplesToInclude)
% This function computes a variance model for each data column (sample) in
% d. The second input samplesToInclude can be used to specify only a subset
% of the data columns that should be used for estimating the average
% frequency in the library of each of the control guides.

if nargin<2
    samplesToInclude=1:size(d.data,2);
end

%% Identify control guides
indSafe=find(strcmp(d.gene,'safe'));
indNone=find(strcmp(d.gene,'none'));
indOtherCont=find(boolRegExp(d.gene,'^CONT_[0-9]+'));
indC=union(indSafe,indNone);
indC=union(indC,indOtherCont);

%% Make a data table for control guides. We will use control guide data to calculate the variance models
cont.data=d.data(indC,:);
cont.sums=sum(d.data(indC,:),2);
cont.means=mean(d.data(indC,:),2);
cont.totReads=sum(cont.data,1);

%% Analyze the control guides over all samples
cont.indepSums=sum(cont.data(:,samplesToInclude),2);
%cont.freqs=cont.indepSums/sum(d.totReads(samplesToInclude));
%cont.adjForTotContEnrich=(cont.totReads./d.totReads)/sum(cont.freqs);  % This is a correction for 
%cont.adjForTotContEnrich=(median(cont.data)./d.totReads)/median(cont.freqs);
cont.freqs=mean(cont.data(:,samplesToInclude)./d.totReads(samplesToInclude),2);
adj1=cont.data./d.totReads./cont.freqs;
effTotReads=nanmedian(adj1).*d.totReads;

expNum=10.^(1.2:0.05:3); %expected read numbers to analyze

expCounts=zeros(size(cont.data));
for j=1:size(cont.data,2)
%    expCounts(:,j)=cont.freqs*d.totReads(j)*cont.adjForTotContEnrich(j);
    expCounts(:,j)=cont.freqs*effTotReads(j);
end
deviations=cont.data-expCounts;
diffs=cell(length(expNum),size(cont.data,2));
for j=1:size(cont.data,2)
    for i=1:length(expNum)
        ind=abs(expCounts(:,j) - expNum(i)) < 0.2*expNum(i);
        diffs{i,j}=deviations(ind,j);
        noiseVals(i,j)=sqrt(mean(diffs{i,j}.^2)/expNum(i));
        varianceVals(i,j)=std(diffs{i,j})^2 * length(samplesToInclude)/(length(samplesToInclude)-1); %sqrt(sum(diffs{i,j}.^2));
        countVals(i,j)=length(diffs{i,j});
    end
end
%%
for j=1:size(d.data,2)
    ind=countVals(:,j)>=25;
    varModels(j)=fitline(log(expNum(ind)),log(varianceVals(ind,j)));
end
