%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% COMPARING FITS+PARAMETER VALUES FOR A SINGLE SESSION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expRef =  '2015-10-30_1_Hopkins';
saveDir = '\\basket.cortexlab.net\homes\peterzh\NeuroGLM\ModelFiles';
g = GLM(expRef);

models = {'ifC','Supersaturation-subset','fullContrasts','fullContrasts-subset','CL+CR','CL+CR-subset','C^N','C^N-subset','C50','C50-subset','C^NL^NR','C^NL^NR-subset'};

%% Fit several models to a single dataset and plot their different fits.
% Note this only works with 1D contrast sessions
for m = 1:length(models)
    g = g.setModel(models{m}).fit([]);
    
    subplot(2,length(models),m);
    g.plotFit;
    subplot(2,length(models),m+length(models));
    g.plotParams;
end

%% Leave-1-out fitting over all the models
g = g.CV; %Set CV partitioning;

for m = 1:length(models)
    g = g.setModel(models{m});
    g = g.fit;
    save(fullfile(saveDir,[expRef '_' g.modelString '_crossvalidated.mat']),'g');
    disp('done');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MODEL COMPARISON OVER MULTIPLE SESSIONS %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expRefs = {'2015-05-28_1_Laveran',...
    '2015-05-29_1_Laveran',...
    '2015-05-30_1_Laveran',...
    '2015-05-31_1_Laveran',...
    '2015-06-01_1_Laveran',...
    '2015-06-02_1_Laveran',...
    '2015-06-03_1_Laveran',...
    '2015-06-04_1_Laveran'};

% expRefs = {'2014-11-30_1_M140528NS1',...
%            '2014-12-01_1_M140528NS1',...
%            '2014-12-02_1_M140528NS1'};
% %
% expRefs = {'2015-01-24_1_M140617NS1',...
%            '2015-01-25_1_M140617NS1'};

% expRefs = {'2015-09-21_1_Hopkins','2015-09-21_2_Eijkman'};

models = {'Offset','ifC','C^N','C^N-subset','C^NL^NR','C^NL^NR-subset','C50','C50-subset','Supersaturation-subset'};

saveDir = '\\basket.cortexlab.net\homes\peterzh\NeuroGLM\ModelFiles';

for b = 1:length(expRefs)
    disp(['File: ' expRefs{b}]);
    %     g = GLM(expRefs{b});
    g = laserGLM(expRefs{b});
    for m = 1:length(models)
        g = g.setModel(models{m});
        g = g.fitCV;
        save(fullfile(saveDir,[expRefs{b} '_' g.modelString '_crossvalidated.mat']),'g');
        disp('done');
    end
end

%% Plot multiple session performances (bits per trial)
saveDir = '\\basket.cortexlab.net\homes\peterzh\NeuroGLM\ModelFiles';

% expRefs = {'2015-05-28_1_Laveran',...
%            '2015-05-29_1_Laveran',...
%            '2015-05-30_1_Laveran',...
%            '2015-05-31_1_Laveran',...
%            '2015-06-01_1_Laveran',...
%            '2015-06-02_1_Laveran',...
%            '2015-06-03_1_Laveran',...
%            '2015-06-04_1_Laveran'};
% %
% expRefs = {'2014-11-30_1_M140528NS1',...
%            '2014-12-01_1_M140528NS1',...
%            '2014-12-02_1_M140528NS1'};
%
% % %
% expRefs = {'2015-09-21_1_Hopkins','2015-09-21_2_Eijkman'};
expRefs = {'2015-10-30_1_Hopkins'};

%
% models = {'Offset','ifC','fullContrasts','fullContrasts-subset','C^N','C^N-subset','C^NL^NR','C^NL^NR-subset','C50','C50-subset','Supersaturation-subset'};
% modelLabels = {'offset','ifC','full','full sub','C^N','C^N sub','C^{NL,NR}','C^{NL,NR} sub','C50','C50 sub','Sup sub'};

models = {'Offset','C^N','C^N-subset','C^NL^NR','C^NL^NR-subset','C50','C50-subset','Supersaturation-subset'};
modelLabels = {'Offset','C^N','C^N sub','C^{NL,NR}','C^{NL,NR} sub','C50','C50 sub','Sup sub'};


BitsperTrial = [];
BaselineEntropy = [];
ModelID = [];
for m = 1:length(models)
    p_hat=[];
    for b = 1:length(expRefs)
        load(fullfile(saveDir,[expRefs{b} '_' models{m} '_crossvalidated.mat']));
        p_hat = g.p_hat;
        nanIdx = find(isnan(p_hat));
        p_hat(nanIdx)=[]; %remove failed predictions
        
        bpt = -sum(log2(p_hat))/length(p_hat);
        BitsperTrial(m,b) = bpt;
        
        %Calculate baseline, if the model was guessing based only on the
        %fact that the mouse has a certain total ratio of L:R:NG
        tab = tabulate(g.data.response);
        tab = tab(:,3)/100;
        BaselineEntropy(m,b) = -sum(tab.*log2(tab));
        
    end
end
figure;
subplot(1,2,1);
bar(BitsperTrial,'grouped','EdgeColor','none');
ylabel('Bits per trial');
set(gca,'XTickLabel',modelLabels);
subplot(1,2,2);
bar(BitsperTrial,'stacked','EdgeColor','none');
ylabel('Sum of (bits per trial) over all sessions');
set(gca,'XTickLabel',modelLabels);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% TESTING WITH SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expRef = '2015-10-30_1_Hopkins'; %Used for getting input contrasts from only
trueModel = 'C^N-subset';
trueModelParams = [-0.3 5.8 -1.3 6.1 0.5];
g = simulateGLM(expRef,trueModel,trueModelParams);
g.parameterFits = g.trueParameters;
g.plotFit;

%% Run CV and save results
saveDir = '\\basket.cortexlab.net\homes\peterzh\NeuroGLM\ModelFiles';
models = {'Offset','ifC','CL+CR','CL+CR-subset','C^N','C^N-subset','C^NL^NR','C^NL^NR-subset','C50','C50-subset','Supersaturation-subset'};
modelLabels = {'Offset','ifC','CL+CR','CL+CR sub','C^N','C^N sub','C^{NL,NR}','C^{NL,NR} sub','C50','C50 sub','Sup sub'};

for m = 1:length(models)
    g = g.setModel(models{m});
    g = g.fitCV(400);
    save(fullfile(saveDir,[g.expRef '_' g.modelString '_crossvalidated.mat']),'g');
    disp('done');
end

%% Plot the crossvalidated results of model fitting
BitsperTrial = [];
BaselineEntropy = [];
ModelID = [];
for m = 1:length(models)
    p_hat=[];
    load(fullfile(saveDir,['Simulation2_' trueModel '_' models{m} '_crossvalidated.mat']));
    p_hat = g.p_hat;
    nanIdx = find(isnan(p_hat));
    p_hat(nanIdx)=[]; %remove failed predictions
    
    bpt = -sum(log2(p_hat))/length(p_hat);
    BitsperTrial(m,1) = bpt;
    
    %Calculate baseline, if the model was guessing based only on the
    %fact that the mouse has a certain total ratio of L:R:NG
    tab = tabulate(g.data.response);
    tab = tab(:,3)/100;
    BaselineEntropy(m,1) = -sum(tab.*log2(tab));
    
end
figure;
bar(BitsperTrial,'grouped','EdgeColor','none');
ylabel('Bits per trial');
set(gca,'XTickLabel',modelLabels);
%ylim([0.9 1]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% COMBINING MULTIPLE SESSIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

