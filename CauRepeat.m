function [timeCau,EEout1, EEout2, coeff, pvalue] = CauRepeat(timeSeriesX,timeSeriesY,order, Nr, Nl)
%CAU Summary of this function goes here
%   Detailed explanation goes here
%timeSeriesX Y are matrix whose every row is a variable.
%the output is the causality Y to X.
%sr is the sampling rate
% EEout1 is the determinant of the error process for X given by the model
% EEout2 is the the final prediction error E given by the model for both [X Y] (the variance
%   estimate of the white noise input to the AR model)

[rX,cX] = size(timeSeriesX);
[rY,cY] = size(timeSeriesY);

[A1,E1] = armorfrepeat(timeSeriesX, Nr, Nl,order);
[A2,E2] = armorfrepeat([timeSeriesX;timeSeriesY], Nr, Nl,order);
EEout = predictionerror(A2, [timeSeriesX;timeSeriesY], Nr, Nl, order);
EEout1 = det(EEout(1:rX, 1:rX));
EEout2 = reshape(E2, 1, (rX+rY)^2);
coeff = reshape(A2, 1, (rX+rY)^2*order);
% [E1,A1] = UCM_input([timeSeriesX;timeSeriesZ],order);
% [E2,A2] = UCM_input([timeSeriesX;timeSeriesY;timeSeriesZ],order);
% R1 = E1(1:rX,1:rX)-E1(1:rX,rX+1:rX)*inv(E1(rX+1:rX,rX+1:rX))*E1(rX+1:rX,1:rX);
% R2 = E2(1:rX,1:rX)-E2(1:rX,rX+rY+1:rX+rY)*inv(E2(rX+rY+1:rX+rY,rX+rY+1:rX+rY))*E2(rX+rY+1:rX+rY,1:rX);
R1 = E1;
R2 = E2(1:rX,1:rX);
timeCau = log(det(R1)/det(R2));

%% for one dimension only !!!!!!
df1 = order;
df2 = cX - order - order - 1;
F = (R1 - R2) / R2 * df2 / df1;
pvalue = 1-fcdf(F,df1,df2);



