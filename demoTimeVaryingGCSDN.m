% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %  Time-varying Granger causality with signal-dependent noise
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Version: 1.0.0.0.1
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Description:
% %        This tool implements the time-varying GCSDN algorithm to 
% %        estimate one causal structure at each local time window in a time-series 
% %        observation of some system by borrowing the strength of the whole time-series. 
% %        Together all the locally estimated causal structure establish a time-varying 
% %        interactions among elements of the system. The interaction between two elements
% %        is measured by both strength and level of significance of the infromation flow 
% %        between these two elements at two directions.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % 
% % Inputs:
% %        RoiData: nVol*nROI matrix, time-series data of ROIs one subject
% %                 (nVol:number of volume,i.e. length of the time-series, nROI: number of ROIs)
% %         IndSeq: [nWin x 2 double], specification of local time windows 
% %                  nWin: number of local time windows; 
% %                  local time windows were defined starting and ending
% %                  indexes in the time-series
% %       extraWei: [nWin x 1 double], user specified weight of each local time windows
% %                 e.g. local time windows with head motion exceeding a given threshold
% %                 can be asigned 0 weights; 
% %                 if the Head motion >threshold, extraWei(k)=1,
% %                 else extraWei(k)=0
% %
% % Output:result
% %        coefEST: cell{nWin,1}, estimated coeff 
% %                 A:coef(1:4,:), C:coef(5:6,:), B:coef(7:10,:)
% %         llRate: [nWin x nEdge double]; likelyhood rate of GCSDN
% %          IFsig: [nWin x nEdge double]; significoncy of information flow
% %             IF: [nWin x nEdge double]; stength of the information flow, GCSDN
% %                 measurementy
% %         stable: [nWin x nEdge double]; 
% %                 -1 mean stsble, 1000 means unstable
% % 
% % Example: 
% %          ROIdata = randn(100,3);
% %           IndSeq = [1 20;21 40;41 60;61 100];
% %         extraWei = [0 0 0 0];
% %           result = demoTimeVaryingGCSDN(ROIdata,IndSeq,extraWei);
% % 
% % Please cite the following reference:
% %   Dynamic effective connectivity reveals an in-default wiring for 
% %   prosociality and an on-demand circuit for deception in the brain, by 
% %   Baobao Pan, Qiang Luo, et al. (submitted), 2018.
% %   
% % Contact:
% %   If you have any question regarding this tool, contact can be
% %   addressed to the following email address:
% %   Qiang Luo, Fudan University, mrqiangluo@gmail.com
% %   The code was last updated by Baobao Pan (pbb_194@163.com) on 16 Mar 2018.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function result = demoTimeVaryingGCSDN(ROIdata,IndSeq,extraWei)

% % - - - - - - - - parameter setting - - - - - - - % %              
nROI = size(ROIdata,2);                % number of ROIs
nWin = size(IndSeq,1);               % number of time windows
ar_order = 1;                          % order of AR model
bekk_order = 1;                        % order of AR-BEKK model
indexX = 1;                            % label of time seriers X
indexY = 2;                            % label of time seriers Y

% % - - - - - for spectral analysis only - - - - - - % %
sr = 0.5;                              % sampling frequency  
fd.EDFreq = sr/2; fd.STFreq = 0;       % parameter for frequency domain
fd.NFFT = 256; fd.fs = sr;

nComb = nchoosek(nROI,2);
nEdge = nComb*2;
combination = nchoosek(1:nROI,2);
timeCauREbekkSig_Sub = nan(nWin,nEdge);
timeCauREbekkSig_tvSub = nan(nWin,nEdge);
timeCauREbekk_tvSub = nan(nWin,nEdge);
stable_tvSub = nan(nWin,nComb);
para = cell(1,1);
trial = reshape(1:nWin,1,nWin)';
nWin = size(trial,1);
exWin = extraWei;
for idTrial = 1:nWin
    for i = 1 :size(combination,1)
        %disp(['Trial: ',num2str(idTrial),'  comb: ',num2str(combination(i,:))]);
        clear input_data
        input_data.timeseriesdata = ROIdata(:,combination(i,:));
        input_data.Nl = IndSeq;
        input_data.Nr = size(input_data.Nl,1);
        input_data.exclude = exWin;
        outputarmabekk = mv_grangerarmabekk4RepeatTimeVary(input_data, trial,...
            idTrial, ar_order, bekk_order, indexX, indexY, fd);
        timeCauREbekkSig_Sub(idTrial, 2*i-1) = outputarmabekk.granger(1);
        timeCauREbekkSig_Sub(idTrial, 2*i) = outputarmabekk.granger(2);
        timeCauREbekkSig_tvSub(idTrial, 2*i-1) = outputarmabekk.granger(3);
        timeCauREbekkSig_tvSub(idTrial, 2*i) = outputarmabekk.granger(4);
        timeCauREbekk_tvSub(idTrial, 2*i-1) = outputarmabekk.granger(5);
        timeCauREbekk_tvSub(idTrial, 2*i) = outputarmabekk.granger(6);
        stable_tvSub(idTrial, i) = stationary_constraint(outputarmabekk.parameters, 1, 1, 2, 1, 1);
        para{idTrial}(:,i) = outputarmabekk.parameters;
    end
end
result.coefEST = para;
result.llRate = timeCauREbekkSig_Sub;
result.IFsig = timeCauREbekkSig_tvSub;
result.IF  =timeCauREbekk_tvSub;
result.stable = stable_tvSub;
