% Y->X
function [timeCau,R1,R2,coeff,pvalue] = Cau(timeSeriesX,timeSeriesY,order)

[rX,cX] = size(timeSeriesX);
[rY,cY] = size(timeSeriesY);

[A1,E1] = armorf(timeSeriesX,1,cX,order);
[A2,E2] = armorf([timeSeriesX;timeSeriesY],1,cX,order);

R1=E1;
R2=E2(1:rX,1:rX);
timeCau = log(trace(R1)/trace(R2));
coeff = A2(1,2);

%% for one dimension only !!!!!!
df1 = order;
df2 = cX - order - order - 1;
F = (R1 - R2) / R2 * df2 / df1;
pvalue = 1-fcdf(F,df1,df2);