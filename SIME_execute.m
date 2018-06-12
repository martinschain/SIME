function[Vnd_Out,Vts_Out,FinalModel,KFinal] = SIME_execute(Data,Settings)

% Executes Simultanous estimation of Vnd for time activity curves stored in
% structure Data, according to the setting specified in structure Settings. 
% The structure Data must contain the following fields:
% Data.name    :        Unique subject name or ID
% Data.ROIData :        Array of regional TACs (frames x ROIs)
% Data.ROInames:        Cell array of ROI names
% Data.MidTime :        Vector of times corresponding to the TACs
% Data.inFcn   :        Vector for the metabolite corrected arterial input function
% Data.wb      :        Vector of radioactivity in whole blood
% Data.inFcn_t :        Vector of times corresponding to inFcn and wb
% Data.Weights :        Time weights to be used when fitting the constrained 2TCM (nbrFramesx1 vextor)
% Data.costFcnWeights:  ROI Weights to be used for calculating the costfunction (nbrROIsx1 vector)
% 
% The structure Settings contains the following fields
% Settings.gridVnd  :     A vector of Vnd values that should be evaluated  
% Settings.vB       :     Fractional blood volume, either a scalar (i.e., 0.05), or, if left empty, vB will be fitted
% Settings.doPlot   :     Set to 1 plot and save output figures
% Settings.pathOut  :     Path to where results are saved.


%% Resample the input function so that it has equal time steps
if Data.inFcn_t(end) > 600
    disp('Guessing input function time is given in seconds, converting to min')
    Data.inFcn_t = Data.inFcn_t/60;
end
stepSize = 1/30; % Sampling distance set to 2 second
t_interp = Data.inFcn_t(1):stepSize:Data.inFcn_t(end);
Data.inFcn = interp1(Data.inFcn_t,Data.inFcn,t_interp);
Data.wb = interp1(Data.inFcn_t,Data.wb,t_interp);
Data.inFcn_t = t_interp;

%% Check time unit for TACs
if Data.MidTime(end) > 600
    disp('Guessing TAC time is given in seconds, converting to min')
    Data.MidTime = Data.MidTime/60;
end

%% For every choice of Vnd, estimate the model
tic
nbrOfROIs = size(Data.ROIData, 2);
Ks_tROIs = zeros(nbrOfROIs,4,length(Settings.gridVnd));
Model2TCC = zeros(size(Data.ROIData,1),size(Data.ROIData,2),length(Settings.gridVnd));
RSStot = zeros(nbrOfROIs,length(Settings.gridVnd));

for gvnd = 1:length(Settings.gridVnd)
    Data.Vnd = Settings.gridVnd(gvnd);
    disp(['Estimating Ks for Vnd = ' num2str(Settings.gridVnd(gvnd)) ', (' num2str(round(100*gvnd/length(Settings.gridVnd))) '% done, ' num2str(toc) 's passed).'])
    for roi = 1:nbrOfROIs
        [modelTAC,K,RSS]      = SIME_convolve(Data.inFcn,Data.wb,Data.inFcn_t,Data.ROIData(:,roi),Data.MidTime,'2TCMfixedVnd',Data.Vnd,Data.Weights,Settings);
        Ks_tROIs(roi,:,gvnd)  = K;
        Model2TCC(:,roi,gvnd) = modelTAC';
        RSStot(roi,gvnd) = RSS;
    end
    Data = rmfield(Data,'Vnd');
end


%% Weight the RSS with cost function weights (one per ROI)
costFcnWeights = repmat(Data.costFcnWeights(:),[1 length(Settings.gridVnd)]);
RSStot_weighted = RSStot.*costFcnWeights;

%% Sum the weighted RSS across ROIs to create cost function
costFcn = sum(RSStot_weighted,1);

%% Find minimum position in costFcn, and create output
[~, costFunMinId] = min(costFcn); %Find the position where the function reaches minimum
Vnd_Out = Settings.gridVnd(costFunMinId);
FinalModel = Model2TCC(:,:,costFunMinId);
KFinal = Ks_tROIs(:,:,costFunMinId);
Vts_Out = Vnd_Out.*(1+KFinal(:,2)./KFinal(:,3));
gridVnd = Settings.gridVnd;
if isfield(Data,'VndTrue')
    VndTrue = Data.VndTrue;
    save([Settings.pathOut filesep Data.name '.mat'],'Vnd_Out','KFinal','FinalModel','Vts_Out','costFcn','gridVnd','VndTrue')
else
    save([Settings.pathOut filesep Data.name '.mat'],'Vnd_Out','KFinal','FinalModel','Vts_Out','costFcn','gridVnd')
end

%% Plot the data
if Settings.doPlot
    %% Plot the fits
    figure('visible','on');
    subplot(311), plot(Data.MidTime,Data.ROIData,'o'), title('All fits'), hold on
    subplot(312), plot(Data.MidTime,Data.ROIData*0,'-'), title('All residuals'), hold on
    for i = 1:size(Model2TCC,3)
        subplot(311), plot(Data.MidTime,Model2TCC(:,:,i),':')
        subplot(312), plot(Data.MidTime,Data.ROIData-Model2TCC(:,:,i),'o')
    end
    
    % Fill in the final model and residual
    subplot(311), plot(Data.MidTime,FinalModel,'-','LineWidth',1.5)
    subplot(312), plot(Data.MidTime,Data.ROIData - FinalModel,'*','MarkerSize',14)
    
    % Plot the cost function
    subplot(313), plot(gridVnd,costFcn,'m')
    saveas(gcf,[Settings.pathOut filesep Data.name '.png'])
end

