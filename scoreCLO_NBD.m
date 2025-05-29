function [e,g] = scoreCLO_NBD(d,comparisons,minCounts,varModels)
% SCORECLO computes scores for elements (guides) and genes from CRISPR screen data

if nargin<3
    comparisons={[1 2]};
end
lenComp=length(comparisons);
if nargin<4 | isempty(minCounts)
    minCounts=10;
end

%Resolution of sampling for computing Bayesian estimates
resolution=20;
effRange=-8:0.01:8;
minSumCountsBayesian=100; %minimum for the mean of the counts per sample to compute the prior distribution for effect size


indSafe=find(strcmp(d.gene,'safe'));
indNone=find(strcmp(d.gene,'none'));
indC=union(indSafe,indNone);

% If not provided, compute the variance models
if nargin<4
    samplesToInclude=unique(cell2mat(comparisons));
    varModels= CLO_compute_variance_models(d,samplesToInclude);
end

%% element-level data
e.gene=d.gene;
e.guideID = d.guideID;

for compNum=1:lenComp
    tic;
    samples=comparisons{compNum};
    if size(samples,2)==2  %pairwise comparison between samples
        freq=sum(d.data(:,samples),2)/sum(d.totReads(samples));
        %effTotReads=(median(d.data(indC,samples))./d.totReads(samples))/median(freq(indC)).*d.totReads(samples);
        %freq=mean(d.data(:,samples)./d.totReads(samples),2);
        effTotReads=nanmedian((d.data(indC,samples)./d.totReads(samples))./freq(indC)).*d.totReads(samples);

        % temp1 will be the probabilities of obtaining the observed results
        % with the null hypothesis that the samples are drawn from the same
        % distribution
        meanVals=max(1,freq*effTotReads); %d.totReads(samples);
        varVals=[exp(varModels(samples(1)).b + varModels(samples(1)).m * log(meanVals(:,1))) exp(varModels(samples(2)).b + varModels(samples(2)).m * log(meanVals(:,2)))];
        varVals(varVals<=meanVals)=1.01*meanVals(varVals<=meanVals);
        temp1=negBinDistPDF(d.data(:,samples),meanVals,varVals);

        % temp2 will be the probabilities for the alternate model (the
        % observations match the true probabilities
        meanVals=max(1,d.data(:,samples));  % Setting in a minimum of 1 avoids NaN and -Inf values
        varVals=[exp(varModels(samples(1)).b + varModels(samples(1)).m * log(meanVals(:,1))) exp(varModels(samples(2)).b + varModels(samples(2)).m * log(meanVals(:,2)))];
        varVals(varVals<=meanVals)=1.01*meanVals(varVals<=meanVals);
        temp2=negBinDistPDF(d.data(:,samples),meanVals,varVals);

        % compute the log odds
        e.lod1(:,compNum)=-1*sum(log10(temp1./temp2),2);
        e.lod1(sum(d.data(:,samples),2)<(minCounts*2),compNum)=nan;  % log-odds with the modeled noise distribution

        %compute empirical p-values based on rank order with the controls
        lodList=e.lod1(indC,compNum);
        lodList(isnan(lodList))=[];
        len=length(lodList);
        e.p1(:,compNum)=arrayfun(@(x) (0.5+sum(lodList>x))/(len+1),e.lod1(:,compNum));
        e.p1(isnan(e.lod1(:,compNum)),compNum)=nan;
        %effSize0=log2(d.data(:,samples(1))./d.data(:,samples(2))) - log2(d.totReads(samples(1))/d.totReads(samples(2)));
        effSize0=log2(d.data(:,samples(1))./d.data(:,samples(2))) - log2(effTotReads(1)/effTotReads(2));
        e.effRaw(:,compNum)=effSize0;
        e.lod2(:,compNum)=e.lod1(:,compNum);
        e.lod2(effSize0<0,compNum)=-1*e.lod2(effSize0<0,compNum);
        lodList=e.lod2(indC,compNum);
        lodList(isnan(lodList))=[];
        len=length(lodList);
        e.p2(:,compNum)=arrayfun(@(x) (0.5+sum(lodList>x))/(len+1),e.lod2(:,compNum));
        e.p2(isnan(e.lod1(:,compNum)),compNum)=nan;

        %generate Bayesian estimate for effect sizes
        [f,xi]=ksdensity(effSize0(sum(d.data(:,samples),2)>minSumCountsBayesian*2),effRange,'Kernel',@(x) tpdf(x,3)); %switched to tdist kernel for tail behavior
        f2=exp(smooth(log(f),0.1,'rloess')); % smoothing the distribution to avoid artifacts from local peaks
        %f2(xi>-1 & xi<1)=f(xi>-1 & xi<1);
        f=f2;
        ratio1=zeros(size(effSize0));
        for i=1:length(effSize0)
            if ~isnan(e.p1(i,compNum)) %don't compute for samples without enough counts
                %p1vals=[d.data(i,samples(1))/d.totReads(samples(1)) freq(i)];
                p1vals=[d.data(i,samples(1))/effTotReads(1) freq(i)];
                if diff(p1vals)~=0
                    p1vals=[p1vals(1):(diff(p1vals)/resolution):p1vals(2)];
                else
                    p1vals=[0.99*p1vals(1) p1vals(1) 1.01*p1vals(1)];
                end
                %p2vals=[d.data(i,samples(2))/d.totReads(samples(2)) freq(i)];
                p2vals=[d.data(i,samples(2))/effTotReads(2) freq(i)];
                if diff(p2vals)~=0
                    p2vals=[p2vals(1):(diff(p2vals)/resolution):p2vals(2)];
                else
                    p2vals=[0.99*p2vals(1) p2vals(1) 1.01*p2vals(1)];
                end
                %likemat=gammaDistPDF(d.data(i,samples(1)),p1vals*d.totReads(samples(1)),theta(samples(1)))'*gammaDistPDF(d.data(i,samples(2)),p2vals*d.totReads(samples(2)),theta(samples(2)));
                meanVals1=p1vals*effTotReads(1); %d.totReads(samples(1));
                meanVals2=p2vals*effTotReads(2); %d.totReads(samples(2));
                varVals1=exp(varModels(samples(1)).b + varModels(samples(1)).m * log(meanVals1));
                varVals1(varVals1<=meanVals1)=1.01*meanVals1(varVals1<=meanVals1);
                varVals2=exp(varModels(samples(2)).b + varModels(samples(2)).m * log(meanVals2));
                varVals2(varVals2<=meanVals2)=1.01*meanVals2(varVals2<=meanVals2);
                likemat=negBinDistPDF(d.data(i,samples(1)),meanVals1,varVals1)'*negBinDistPDF(d.data(i,samples(2)),meanVals2,varVals2);
                ratiomat=log2(p1vals'*(1./p2vals));
                priormat=interp1(xi,f,ratiomat);
                likematPost=likemat.*priormat;
                [i1,i2]=find(likematPost==max(likematPost(:)),1);
                ratio1(i)=ratiomat(i1,i2);
            end
        end
        e.effNorm(:,compNum)=ratio1;
    else % 4 sample comparison, such as for endocytosis screen
        freq=sum(d.data(:,samples),2)./sum(d.totReads(samples));
        effTotReads=nanmedian((d.data(indC,samples)./d.totReads(samples))./freq(indC)).*d.totReads(samples);

        % pH=sum(d.data(:,samples([1 3])),2)/sum(d.totReads(samples([1 3]))); % read fraction in the high bin samples
        % pL=sum(d.data(:,samples([2 4])),2)/sum(d.totReads(samples([2 4]))); % read fraction in the low bin samples
        % pR=sqrt((sum(d.data(:,samples([1 2])),2)/sum(d.totReads(samples([1 2]))))./(sum(d.data(:,samples([3 4])),2)/sum(d.totReads(samples([3 4]))))); %corrects for unequal total read frequencies for each guide for the stimulated vs unstimulated screens
        pH=sum(d.data(:,samples([1 3])),2)/sum(effTotReads(samples([1 3]))); % read fraction in the high bin samples
        pL=sum(d.data(:,samples([2 4])),2)/sum(effTotReads(samples([2 4]))); % read fraction in the low bin samples
        pR=sqrt((sum(d.data(:,samples([1 2])),2)/sum(effTotReads(samples([1 2]))))./(sum(d.data(:,samples([3 4])),2)/sum(effTotReads(samples([3 4]))))); %corrects for unequal total read frequencies for each guide for the stimulated vs unstimulated screens
        indFit=min(d.data(:,samples),[],2)>minCounts;
        
        % Smooth model for systematic trends in read frequencies
        %meansOld=[pH.*pR pL.*pR pH./pR pL./pR].*repmat(d.totReads(samples),[size(pH,1) 1]);  % Our older method of estimation, without any trend correction
        meansOld=[pH.*pR pL.*pR pH./pR pL./pR].*repmat(effTotReads(samples),[size(pH,1) 1]);  % Our older method of estimation, without any trend correction
        for j=1:4         
            [bins1,y1]=evalBySlidingBinsFixedN(log(pH(indFit)./pL(indFit)),log((d.data(indFit,samples(j))/effTotReads(j))./(meansOld(indFit,j)/effTotReads(j))),200,'median',50);
            y2=smooth(bins1,y1,0.2,'loess');
            [~,indUnique]=unique(bins1);
            bins1=bins1(indUnique);
            y2=y2(indUnique);
            trendAdj=interp1(bins1,y2,log(pH./pL),'linear',NaN);
            trendAdj(log(pH./pL)<min(bins1))=y2(1);
            trendAdj(log(pH./pL)>max(bins1))=y2(end);         % fill in data for extreme points using the correction for the closest data point   
            trendAdj=exp(trendAdj);
            meanVals(:,j)=meansOld(:,j).*trendAdj;
        end


        % temp1 will be the null model (same ratio in high vs low
        % fractions, after correcting for overall abundance/frequency in
        % screen1 vs screen2 (stimulated vs unstimulated)
        %meanVals=[pH.*pR pL.*pR pH./pR pL./pR].*repmat(d.totReads(samples),[size(pH,1) 1]);
        %meanVals=[pH.*pR./pT pL.*pR.*pT pH./pR.*pT pL./pR./pT].*repmat(d.totReads(samples),[size(pH,1) 1]);
        varVals=exp(repmat([varModels(samples).b],size(pH)) + repmat([varModels(samples).m],size(pH)).*log(meanVals));
        varVals(varVals<=meanVals)=1.01*meanVals(varVals<=meanVals);
        temp1=negBinDistPDF(d.data(:,samples),meanVals,varVals);
        
        % temp2 will be the alternate model
        meanVals=max(1,d.data(:,samples));   % Setting in a minimum of 1 avoids NaN and -Inf values
        varVals=exp(repmat([varModels(samples).b],size(pH)) + repmat([varModels(samples).m],size(pH)).*log(meanVals));
        varVals(varVals<=meanVals)=1.01*meanVals(varVals<=meanVals);
        temp2=negBinDistPDF(d.data(:,samples),meanVals,varVals);
        
        e.lod1(:,compNum)=-1*sum(log10(temp1./temp2),2);
        e.lod1(sum(d.data(:,samples),2)<(minCounts*length(samples)),compNum)=nan;  % log-odds with the modeled noise distribution

        %compute empirical p-values based on rank order with the controls
        lodList=e.lod1(indC,compNum);
        lodList(isnan(lodList))=[];
        len=length(lodList);
        e.p1(:,compNum)=arrayfun(@(x) (0.5+sum(lodList>x))/(len+1),e.lod1(:,compNum));
        e.p1(isnan(e.lod1(:,compNum)),compNum)=nan;
        effSize0=log2(d.data(:,samples(1))./d.data(:,samples(2))./d.data(:,samples(3)).*d.data(:,samples(4))) - log2(d.totReads(samples(1))/d.totReads(samples(2))/d.totReads(samples(3))*d.totReads(samples(4)));
        e.effRaw(:,compNum)=effSize0;
        e.lod2(:,compNum)=e.lod1(:,compNum);
        e.lod2(effSize0<0,compNum)=-1*e.lod2(effSize0<0,compNum);
        lodList=e.lod2(indC,compNum);
        lodList(isnan(lodList))=[];
        len=length(lodList);
        e.p2(:,compNum)=arrayfun(@(x) (0.5+sum(lodList>x))/(len+1),e.lod2(:,compNum));
        e.p2(isnan(e.lod1(:,compNum)),compNum)=nan;

        %generate Bayesian estimate for effect sizes
        [f,xi]=ksdensity(effSize0(sum(d.data(:,samples),2)>minSumCountsBayesian*4),effRange,'Kernel',@(x) tpdf(x,3)); %switched to tdist kernel for tail behavior
        f2=exp(smooth(log(f),0.1,'rloess')); % smoothing the distribution to avoid artifacts from local peaks
        %f2(xi>-1 & xi<1)=f(xi>-1 & xi<1);
        f=f2;
        ratio1=zeros(size(effSize0));
        for i=1:length(effSize0)
            if ~isnan(e.p1(i,compNum)) %don't compute for samples without enough counts
                p1vals=[d.data(i,samples(1))/d.totReads(samples(1)) pH(i)*pR(i)];
                p1vals=[p1vals(1):(diff(p1vals)/resolution):p1vals(2)];
                p2vals=[d.data(i,samples(2))/d.totReads(samples(2)) pL(i)*pR(i)];
                p2vals=[p2vals(1):(diff(p2vals)/resolution):p2vals(2)];
                p3vals=[d.data(i,samples(3))/d.totReads(samples(3)) pH(i)/pR(i)];
                p3vals=[p3vals(1):(diff(p3vals)/resolution):p3vals(2)];
                p4vals=[d.data(i,samples(4))/d.totReads(samples(4)) pL(i)/pR(i)];
                p4vals=[p4vals(1):(diff(p4vals)/resolution):p4vals(2)];
                meanVals1=p1vals*d.totReads(samples(1));
                meanVals2=p2vals*d.totReads(samples(2));
                meanVals3=p3vals*d.totReads(samples(3));
                meanVals4=p4vals*d.totReads(samples(4));
                varVals1=exp(varModels(samples(1)).b + varModels(samples(1)).m * log(meanVals1));
                varVals1(varVals1<=meanVals1)=1.01*meanVals1(varVals1<=meanVals1);
                varVals2=exp(varModels(samples(2)).b + varModels(samples(2)).m * log(meanVals2));
                varVals2(varVals2<=meanVals2)=1.01*meanVals2(varVals2<=meanVals2);
                varVals3=exp(varModels(samples(3)).b + varModels(samples(3)).m * log(meanVals3));
                varVals3(varVals3<=meanVals3)=1.01*meanVals3(varVals3<=meanVals3);
                varVals4=exp(varModels(samples(4)).b + varModels(samples(4)).m * log(meanVals4));
                varVals4(varVals4<=meanVals4)=1.01*meanVals4(varVals4<=meanVals4);
                likemat12=negBinDistPDF(d.data(i,samples(1)),meanVals1,varVals1)'*negBinDistPDF(d.data(i,samples(2)),meanVals2,varVals2);
                likemat34=negBinDistPDF(d.data(i,samples(3)),meanVals3,varVals3)'*negBinDistPDF(d.data(i,samples(4)),meanVals4,varVals4);
                likemat=outerProduct(likemat12,likemat34);
                ratiomat=log2(outerProduct(p1vals'*(1./p2vals),(1./p3vals)'*p4vals));
                priormat=interp1(xi,f,ratiomat);
                likematPost=likemat.*priormat;
                [i1,i2]=find(likematPost==max(likematPost(:)),1);
                ratio1(i)=ratiomat(i1,i2);
            end
        end
        e.effNorm(:,compNum)=ratio1;
    end
    toc
    fprintf('Finished element-level scoring for comparison %i\n',compNum);
end
% Apply max/min bounds and get rid of any infinite values
e.lod2(e.lod2<-100)=-100;
e.lod2(e.lod2>100)=100;

% Gene-level data
g.ensg=setdiff(d.ensg,{'0None','0Safe'});
[~,ia,ib]=intersect(d.ensg,g.ensg);
g.gene=d.gene(ia);
%g.gene=unique(d.gene);
g.p2=nan(length(g.gene),lenComp);
g.p2rev=nan(size(g.p2));
g.lod2=nan(size(g.p2));
g.eff=nan(size(g.p2));
for compNum=1:lenComp
    tic
    for i=1:length(g.ensg)
        ind=strcmp(g.ensg(i),d.ensg) & ~isnan(e.lod1(:,compNum));
        if sum(ind)>0
            g.p2(i,compNum)=combinepvalues(e.p2(ind,compNum));
            g.p2rev(i,compNum)=combinepvalues(1-e.p2(ind,compNum));
            g.lod2(i,compNum)=sum(e.lod2(ind,compNum));
            if g.lod2(i,compNum)>0
                g.eff(i,compNum)=max(e.effNorm(ind,compNum));
            else
                g.eff(i,compNum)=min(e.effNorm(ind,compNum));
            end
        end
    end
    toc;
    fprintf('Finished gene-level scoring for comparison %i\n',compNum);
end
g.lp2=-1*log10(g.p2);  %positive sign for low p-values
g.lp2(g.p2>g.p2rev)=log10(g.p2rev(g.p2>g.p2rev));  %negative sign for high p-values (close to 1)
g.eff2=g.eff.*max(0,(1-2*min(g.p2,g.p2rev))); % shrink effect size estimates for genes with high p-values
