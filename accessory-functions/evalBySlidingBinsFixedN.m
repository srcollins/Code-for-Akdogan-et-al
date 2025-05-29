function [xOut,yOut]=evalBySlidingBinsFixedN(xIn,yIn,n,func,step)
%Evaluates the function func on a sliding bin of width n (n being a number
%of data points) over the paired data (xIn,yIn)

if nargin<5
    step=1;
end
temp=[xIn(:) yIn(:)];
temp(isnan(xIn(:)+yIn(:)),:)=[];
temp=sortrows(temp,1);
xIn=temp(:,1); yIn=temp(:,2);
for i=1:((length(xIn)-n)/step +1);
    start=((i-1)*step+1);
    range=start:(start+n-1);
    xOut(i)=median(xIn(range));
    yOut(i)=feval(func,yIn(range));
end
