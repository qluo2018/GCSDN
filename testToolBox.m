clear;close all;clc;

% % A simple example data 
% ROIdata = randn(100,3);
% IndSeq = [1 40;41 100];
% extraWei = [0 0];
% result = demoTimeVaryingGCSDN(ROIdata,IndSeq,extraWei);
% disp('The End ... ...');


% % % - - - -example for parallel computing
load('exampleData.mat');    % Real BOLD data of 3 subjects
nSub = length(postROI);
result = cell(nSub,1);
parfor subj = 1:nSub  
   ROIdata = postROI{subj}; 
   IndSeq = postInd{subj};
   extraWei = exclude(subj,:);
   result{subj,1} = demoTimeVaryingGCSDN(ROIdata,IndSeq,extraWei);
   disp(['subject: ',num2str(subj)])
end
