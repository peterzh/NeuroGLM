classdef neurGLM
    %For analysis of choiceworld neural data using the GLM
    
    properties (Access=public)
        name;
        epochs = {'stimulusCueStarted',[-0.5 1.5];
            'stimulusCueEnded',[-0.5 1.5];
            'goCue',[-1 1];
            'responseMoveStartTimes',[-1 +1];
            'responseMade',[-1.5 0.5];
            'feedbackStarted',[-1 1]};
        cwLabels;
        cwEvents;
        spikes;
        timebinWidth=0.01;
        split;
    end
    
    properties (Access=private)
        numEval=100;
        numCVfolds=5;
    end
    
    methods
        function obj = neurGLM(cw_dir,name)
            obj.name = name;
            disp('Loading data...');
            cwl=load(fullfile(cw_dir,'cwLabels.mat'));
            s=load(fullfile(cw_dir,'spikes.mat'));
            cwe=load(fullfile(cw_dir,'cwEvents.mat'));
            inclTrials = [cwl.cwLabels.repeatNum==1 & cwl.cwLabels.inclTrials]';
%             inclTrials = ones(length(cwl.cwLabels.repeatNum),1);
            
            %remove invalid trials
            cwl.cwLabels=rmfield(cwl.cwLabels,'contrastIndVals');
            obj.cwLabels=structfun(@(a)(a(inclTrials)'),cwl.cwLabels,'uni',0);
            obj.cwEvents=structfun(@(a)(a(inclTrials)),cwe.cwEvents,'uni',0);
            obj.spikes=s.spikes;
            
            
            obj.split.Groups = {obj.spikes.shanks.singleUnitInds'};
            obj.split.Type = 'none';
            obj.split.splitVals = {''};
            
            if any(obj.cwLabels.responseMade==3)
                obj = obj.NoGoResampling;
            end
        end
        
        function obj = calculatePopRate(obj,varargin)
            disp('Calculating population rate...');
            
            obj.spikes.shanks.populationRateT = 0:obj.timebinWidth:obj.cwEvents.trialEnded(end);
           
            
            numUnits = cellfun(@length,obj.split.Groups);
            unitIDs = obj.split.Groups;
            obj.spikes.shanks.populationRate = zeros(length(obj.split.Groups),length(obj.spikes.shanks.populationRateT));
            
            for s = 1:length(unitIDs)
                sTimes_all=arrayfun(@(u){obj.spikes.shanks.units(u).spikeTimes},unitIDs{s});
                sTimes_all(cellfun(@isempty,sTimes_all))=[];
                sTimes=cell2mat(sTimes_all);
                
                sTimes(sTimes>obj.cwEvents.trialEnded(end))=[];
                popRate = hist(sTimes,obj.spikes.shanks.populationRateT);
                obj.spikes.shanks.populationRate(s,:) = popRate;
                
            end
            
            figure;
            plot(obj.spikes.shanks.populationRateT',obj.spikes.shanks.populationRate');
            xlabel('Time (sec)');
            ylabel('PopRate');
            if ~isempty(obj.split.Groups)
                legend(obj.split.splitVals)
            end
            
            if ~isempty(obj.split.Groups)
                filename = [obj.name '_PopRateGroups'];
                
                set(get(gcf,'children'),'fontsize',15);
                figDir = 'B:\figures\GLM+NeuralActivity';
                savefig(fullfile(figDir,[filename '.fig']));
                print(fullfile(figDir,[filename '.pdf' ]),'-dpdf','-painters');
            end
            
        end
        
        function obj = setSplit(obj,what)
            %returns a neurGLM object containing only spiking data within a
            %depth range
            switch(what)
                case 'none'
                    obj.split=struct;
                    obj.split.Groups = {obj.spikes.shanks.singleUnitInds'};
                    obj.split.Type = 'none';
                    obj.split.splitVals = {''};
                    obj = obj.calculatePopRate;
                case 'ypos'
                    sortedVariable = [obj.spikes.shanks.units(obj.spikes.shanks.singleUnitInds).wfYPosition]';
                    cellXpos = [obj.spikes.shanks.units(obj.spikes.shanks.singleUnitInds).wfXPosition]';
                    meanFR = [obj.spikes.shanks.units(obj.spikes.shanks.singleUnitInds).meanFR]';
                    f=figure;
                    scatter(cellXpos,sortedVariable,[],meanFR,'filled');
                    ylabel('Depth');
%                     set(gca,'Ydir','reverse')
                    [sortedVariable,sortOrder] = sort(sortedVariable);
                    sortedCellIds = obj.spikes.shanks.singleUnitInds(sortOrder);
                    title('Select divisions, RIGHT CLICK to end');
                    hold on;
                    exitFlag = 0;
                    Ydivs = [];
                    while exitFlag == 0
                        [~,y,button] = ginput(1);
                        if button == 3 || length(Ydivs) > 10
                            exitFlag=1;
                        else
                            Ydivs = [Ydivs;y];
                            line([min(cellXpos) max(cellXpos)],[y y]);
                        end
                    end
                    Ydivs = sort(Ydivs);
                    close(f);
                    
                    Ydivs = round(Ydivs,0);
                    obj.split.Groups={};
                    obj.split.splitVals={};
                    for i = 1:length(Ydivs)
                        lidx=sortedVariable <= Ydivs(i);
                        obj.split.Groups{i}=sortedCellIds(lidx)';
                        sortedCellIds(lidx)=nan;
                        sortedVariable(lidx)=nan;
                        
                        if i == 1
                            obj.split.splitVals{i} = [ '0 < ' what ' < ' num2str(Ydivs(i)) ];
                        else
                            obj.split.splitVals{i} = [ num2str(Ydivs(i-1)) ' < ' what ' < ' num2str(Ydivs(i)) ];
                        end
                    end
                    obj.split.Groups{i+1} = sortedCellIds(~isnan(sortedCellIds))';
                    obj.split.splitVals{i+1} = [ what ' > ' num2str(Ydivs(i)) ];
                    
                    %wipe the PR since they now need to be recalculated
                    obj.spikes.shanks.populationRate=[];
                    obj.spikes.shanks.populationRateT=[];
                    obj.split.Type = what;
                    
                    obj = obj.calculatePopRate;
                    
                case 'meanFR'
                    sortedVariable = [obj.spikes.shanks.units(obj.spikes.shanks.singleUnitInds).meanFR]';
                    [sortedVariable,sortOrder] = sort(sortedVariable);
                    sortedCellIds = obj.spikes.shanks.singleUnitInds(sortOrder);
                    
                    f=figure;
                    scatter(ones(length(sortedVariable),1),sortedVariable,'filled');
                    ylabel('Mean FR');
                    title('Select divisions, RIGHT CLICK to end');
                    hold on;
                    exitFlag = 0;
                    Ydivs = [];
                    while exitFlag == 0
                        [~,y,button] = ginput(1);
                        if button == 3 || length(Ydivs) > 10
                            exitFlag=1;
                        else
                            Ydivs = [Ydivs;y];
                            line([0 2],[y y]);
                        end
                    end
                    Ydivs = sort(Ydivs);
                    close(f);
                    
                    Ydivs = round(Ydivs,0);
                    obj.split.Groups={};
                    obj.split.splitVals={};
                    for i = 1:length(Ydivs)
                        lidx=sortedVariable <= Ydivs(i);
                        obj.split.Groups{i}=sortedCellIds(lidx)';
                        sortedCellIds(lidx)=nan;
                        sortedVariable(lidx)=nan;
                        
                        if i == 1
                            obj.split.splitVals{i} = [ '0 < ' what ' < ' num2str(Ydivs(i)) ];
                        else
                            obj.split.splitVals{i} = [ num2str(Ydivs(i-1)) ' < ' what ' < ' num2str(Ydivs(i)) ];
                        end
                    end
                    obj.split.Groups{i+1} = sortedCellIds(~isnan(sortedCellIds))';
                    obj.split.splitVals{i+1} = [ what ' > ' num2str(Ydivs(i)) ];
                    
                    %wipe the PR since they now need to be recalculated
                    obj.spikes.shanks.populationRate=[];
                    obj.spikes.shanks.populationRateT=[];
                    obj.split.Type = what;
                    
                    obj = obj.calculatePopRate;
                    
                otherwise
                    error('No valid split chosen');
            end
            
            
            
        end
        
        function [psthSm,binsSm] = calculatePSTH(obj,neuronID,epochID,splitBy,smoothing)
            epoch = obj.epochs{epochID,1};
            time = obj.epochs{epochID,2};
            
            if strcmp(neuronID,'all')
                neuronID = obj.spikes.shanks.singleUnitInds;
            end
            
            if size(neuronID,2)>size(neuronID,1)
                neuronID = neuronID';
            end
            
            sTimes_all=arrayfun(@(u){obj.spikes.shanks.units(u).spikeTimes},neuronID);
            sTimes_all(cellfun(@isempty,sTimes_all))=[];
            spikeTimes=sort(cell2mat(sTimes_all));
            
            if isempty(spikeTimes)
                spikeTimes = 0;
            end
            
            switch(splitBy)
                case 'choice' %{L,R} or {L,R,NG}
                    trGroups = obj.cwLabels.responseMade;
                case 'stimulus' %{L,R,0}
                    c = obj.cwLabels.contrastRight - obj.cwLabels.contrastLeft;
                    trGroups = nan(length(c),1);
                    trGroups(c==0)=3;
                    trGroups(c<0)=1;
                    trGroups(c>0)=2;
                case 'stim+choice'
                    c = obj.cwLabels.contrastRight - obj.cwLabels.contrastLeft;
                    stim = nan(length(c),1);
                    stim(c==0)=3;
                    stim(c<0)=1;
                    stim(c>0)=2;
                    
                    choice = obj.cwLabels.responseMade;
                    
                    trGroups = [choice stim];
            end
            
            eventTimes = obj.cwEvents.(epoch);
            
            nRanges = size(eventTimes,1);
            ranges = [eventTimes+time(1) eventTimes+time(2)];
            rangeLabel = 1:nRanges;
            
            [groupings,~,ic] = unique(trGroups,'rows');
            numGroupings = size(groupings,1);
            rasterX = cell(1,numGroupings);
            rasterY = cell(1,numGroupings);
            for tr = 1:numGroupings
                %     eventTimes_tr = eventTimes(trGroups==groupings(tr));
%                 ranges_tr = ranges(trGroups==groupings(tr),:);
%                 rangeLabel_tr = rangeLabel(trGroups==groupings(tr));
                ranges_tr = ranges(ic==tr,:);
                
                if size(ranges_tr,1) < 10
                    warning(['There are only ' num2str(size(ranges_tr,1)) ' trials of condition: '  num2str(tr)]);
                end
                
                rangeLabel_tr = rangeLabel(ic==tr);
                
                wr_tr = WithinRanges(spikeTimes, ranges_tr, rangeLabel_tr', 'vector');
                
                stIn = spikeTimes(wr_tr>0); wr_tr = wr_tr(wr_tr>0);
                stRelToEvent = stIn-eventTimes(wr_tr); % subtract the event time corresponding to each spike
                
                [psth(tr,:),bins(:,tr)] = hist(stRelToEvent, time(1):obj.timebinWidth:time(2));
                [rasterX{tr},yy] = rasterize(stRelToEvent);
                rasterY{tr} = yy+reshape(repmat(wr_tr,3,1),1,length(wr_tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
                rasterX{tr} = rasterX{tr}(1:3:end);
                rasterY{tr} = rasterY{tr}(1:3:end);

                %normalise each PSTH by the number of trials used to
                %generate that PSTH
                psth(tr,:) = psth(tr,:)/(size(ranges_tr,1)*obj.timebinWidth);
                
                
                %smoothing
                % smooth the PSTH
                if smoothing>0
                    numEle = sum(smoothing > cumsum(diff(bins(:,tr))));
                    gw = gausswin(numEle,3);
                    smWin = gw./sum(gw);
                    psthSm(tr,:) = conv(psth(tr,:), smWin, 'valid');
                    
                    trimStart = round(numEle/2);
                    trimEnd = numEle - trimStart;
                    binsSm(:,tr) = bins(trimStart : (end-trimEnd),tr);
                else
                    psthSm(tr,:) = psth(tr,:);
                    binsSm(:,tr) = bins(:,tr);
                end
                
            end
            
%             scatter(rasterX{1},rasterY{1},'.')
        end
               
        function [psthSm,binsSm] = PSTH(obj,neuronID,epochID,splitBy,smoothing,legendFlag,errorFlag)
            [psthSm,binsSm] = obj.calculatePSTH(neuronID,epochID,splitBy,smoothing);
            
            if errorFlag == 1
                psthError = obj.calculatePSTHError(neuronID,epochID,splitBy,smoothing);
            end
            
            switch(splitBy)
                case 'choice'
                    plot(binsSm, psthSm');
                    
                    if legendFlag==1
                        if max(obj.cwLabels.responseMade)==2
                            legend('left choice','right choice');
                        elseif max(obj.cwLabels.responseMade)==3
                            legend('left choice','right choice','nogo choice');
                        end
                    end
                case 'stimulus'
                    plot(binsSm, psthSm');
                    if legendFlag==1
                        legend('left contrast','right contrast','zero contrast');
                    end
                case 'stim+choice'

                    if max(obj.cwLabels.responseMade)==2
                        plot(binsSm(:,1),psthSm(1,:)','b-',...
                            binsSm(:,2),psthSm(2,:)','r-',...
                            binsSm(:,3),psthSm(3,:)','k-',...
                            binsSm(:,4),psthSm(4,:)','b--',...
                            binsSm(:,5),psthSm(5,:)','r--',...
                            binsSm(:,6),psthSm(6,:)','k--','LineWidth',0.5);
                        set(gca,'box','off');
                        if legendFlag==1
                            legend('left choice + left contrast',...
                                'left choice + right contrast',...
                                'left choice + zero contrast',...
                                'right choice + left contrast',...
                                'right choice + right contrast',...
                                'right choice + zero contrast',...
                                'Orientation','horizontal');
                        end
                    elseif max(obj.cwLabels.responseMade)==3
                        plot(binsSm(:,1),psthSm(1,:)','b-',...
                            binsSm(:,2),psthSm(2,:)','r-',...
                            binsSm(:,3),psthSm(3,:)','k-',...
                            binsSm(:,4),psthSm(4,:)','b--',...
                            binsSm(:,5),psthSm(5,:)','r--',...
                            binsSm(:,6),psthSm(6,:)','k--',...
                            binsSm(:,7),psthSm(7,:)','b-.',...
                            binsSm(:,8),psthSm(8,:)','r-.',...
                            binsSm(:,9),psthSm(9,:)','k-.','LineWidth',0.5);
                        set(gca,'box','off');
                        if legendFlag==1
                            legend('left choice + left contrast',...
                                'left choice + right contrast',...
                                'left choice + zero contrast',...
                                'right choice + left contrast',...
                                'right choice + right contrast',...
                                'right choice + zero contrast',...
                                'nogo choice + left contrast',...
                                'nogo choice + right contrast',...
                                'nogo choice + zero contrast',...
                                'Orientation','vertical');
                        end
                    end
            end
            xlim(obj.epochs{epochID,2});
            set(gcf,'color','w');
        end
        
        function PSTHError(obj,neuronID,epochID,splitBy,smoothing)
            epoch = obj.epochs{epochID,1};
            time = obj.epochs{epochID,2};
            
            if strcmp(neuronID,'all')
                neuronID = obj.spikes.shanks.singleUnitInds;
            end
            
            if size(neuronID,2)>size(neuronID,1)
                neuronID = neuronID';
            end
            
%             psths = nan(length(neuronID),length(time(1):0.01:time(2)));
            for n = 1:length(neuronID)
                [psths(:,:,n),bins] = obj.calculatePSTH(neuronID(n),epochID,splitBy,smoothing);
            end
            psthMean = sum(psths,3);
            psthError = std(psths,[],3);
            
            
            errH = psthMean+psthError;
            errL = psthMean-psthError;
                    
            figure;
            ax = subplot(1,1,1);
            cols = {'r','g','b','k','m','c'};
            hold on;
            for i = 1:size(psthMean,1)
%             for i = [2,4,6]
                plotWithErrUL(ax,bins(:,i)',psthMean(i,:),[errH(i,:); errL(i,:)],cols{i});
            end
            hold off;
                    
        end
        
        function popPSTH(obj,splitBy,smoothing)
            epo = 1:size(obj.epochs,1);
            
            numgroups = length(obj.split.Groups);
            units = obj.split.Groups;
            
            figure;
            
            a=1;
            for g = 1:numgroups
                cellIdxs = units{g};
%                 figure('name',['Group ' num2str(g) '  ' obj.split.splitVals{g} '  n=' num2str(length(cellIdxs))]); 
                
                for i = 1:length(epo)
                    
                    ax(i)=subplot(numgroups,length(epo),a);
                    
                    
                    obj.PSTH(cellIdxs,i,splitBy,smoothing,(i==length(epo) && g == numgroups),0);

                    if g == 1
                        title([obj.epochs{i,1}]);
                    end
                    
                    if i == 1
                        ylabel(sprintf(['Group ' num2str(g) '  \n' obj.split.splitVals{g} ' \n n=' num2str(length(cellIdxs))]),'FontSize',8);
                    end
                    a=a+1;
                end
                
                linkaxes(ax,'y');
                ylim('auto');
            end
            
            set(gcf,'color','w');
        end
        
        function neurPSTH(obj,neuronIDs,epochID,splitBy,smoothing)
            if strcmp(neuronIDs,'all')
                neuronIDs = obj.spikes.shanks.singleUnitInds;
            end
            numN = length(neuronIDs);
            
            %reorder cells by depth
            try
                cellYpos = [obj.spikes.shanks.units(neuronIDs).wfYPosition]';
                if length(cellYpos) < numN
                    error('Some specified cells did not have a position');
                end
                [~,sortOrder] = sort(cellYpos);
                neuronIDs = neuronIDs(sortOrder);
            catch
                warning('problem reordering by depth')
            end
            
%             psthImageT = [obj.epochs{epochID,2}(1):0.01:obj.epochs{epochID,2}(2)];
%             psthImage = nan(numN,length(psthImageT),1);
            figure;
            for i = 1:numN
                subplot(ceil(sqrt(numN)),ceil(sqrt(numN)),i);
                [psthSm,binsSm]=obj.PSTH(neuronIDs(i),epochID,splitBy,smoothing,0,0);
%                 psthImageDiff(i,:,1) = (psthSm(1,:)-psthSm(4,:))/max(max(psthSm));
                %                 psthImageDiff(i,:,2) = (psthSm(1,:)-psthSm(3,:))/max(max(psthSm));
                %                 psthImageDiff(i,:,3) = (psthSm(2,:)-psthSm(3,:))/max(max(psthSm));
                psthImage(i,:) = psthSm(1,:)/max(max(psthSm));
                
                set(gca,'YTickLabel',{},'box','off','YColor',[1 1 1])
                line([0 0],get(gca,'YLim'),'Color',[0 0 0 0.2]);
                set(gca,'XTickLabel',{},'XColor',[1 1 1])
                
            end
            %             ylim('auto');
            %             linkaxes(ax,'y');
            
            set(gcf,'NextPlot','add');
            axes;
            h = title(['Epoch: ' obj.epochs{epochID,1} ' (' num2str(obj.epochs{epochID,2}) '). SplitBy: ' splitBy]);
            set(gca,'Visible','off');
            set(h,'Visible','on');
            set(gcf,'color','w');
            
            figure;
            
%             psthImageDiff(isnan(psthImageDiff))=0;
            %             for j = 1:3
            %                 subplot(1,3,j);
            imagesc(binsSm(:,1),neuronIDs,psthImage);
            
            
        end
        
        function GPFA_TODO(obj)
            %Gaussian process factor analysis using toolbox from Yu
            %https://users.ece.cmu.edu/~byronyu/software.shtml
        end
        
        function PSTHandLogLik_TODO(obj,epochID,varargin)
            epoch = obj.epochs{epochID,1};
            time = obj.epochs{epochID,2};
            num_choices = max(obj.cwLabels.responseMade);
            
            
            if isempty(obj.split.Groups)
                numGroups = 1;
            else
                numGroups = length(obj.split.Groups);
            end
            
            D = struct;
            D.contrast_cond = [obj.cwLabels.contrastLeft obj.cwLabels.contrastRight];
            D.response = obj.cwLabels.responseMade';
            D.repeatNum = obj.cwLabels.repeatNum';
            D.feedbackType = obj.cwLabels.feedbackType';
            D.eventTimes = obj.cwEvents.(epoch);
            D = getrow(D,obj.inclTrials);
            
            timeseries = linspace(time(1),time(2),obj.numEval);
            popRate = nan(num_choices,length(timeseries));
            popRate_err = nan(num_choices,length(timeseries));
            
            try
                PR = obj.spikes.shanks.populationRate(splitFlag,:);
                PRt = obj.spikes.shanks.populationRateT;
            catch
                reply = input('Population rate not found. Calculate? Y/N [Y]:','s');
                if isempty(reply) || strcmp(reply,'Y') || strcmp(reply,'y')
                    obj = obj.calculatePopRate;
                    PR = obj.spikes.shanks.populationRate(splitFlag,:);
                    PRt = obj.spikes.shanks.populationRateT;
                end
            end
            
            %             f=figure;
            if isempty(obj.split.Groups)
                Units = obj.spikes.shanks(1).singleUnitInds;
            else
                Units = obj.split.Groups{splitFlag};
            end
            numUnits = length(Units);
            
            D.firingRate = nan(length(D.response),numUnits); %firing rates for indiv neurons
            LL_PR=nan(1,length(timeseries)); %loglik for popRate model
            LL_NR=nan(1,length(timeseries)); %loglik for allneur model
            LL_baseline=nan(1,length(timeseries)); %loglik for baseline model
            
            %Calculate delays from all other epochs to overlay..
            all = 1:size(obj.epochs,1);
            other_epochs = all(all~=epochID);
            delays = cell(length(other_epochs),1);
            for i = 1:length(other_epochs)
                epochquery = obj.epochs{other_epochs(i),1};
                delays{i} = obj.cwEvents.(epochquery) - obj.cwEvents.(epoch);
            end
            
            for t = 1:length(timeseries)
                %                 disp([num2str(t) '/' num2str(obj.numEval)]);
                
                %Extract the population rate at this shifted timepoint for all trials
                shiftTime = timeseries(t);
                D.popRates = arrayfun(@(id)PR(find(PRt>=(D.eventTimes(id)+shiftTime),1,'first')),1:length(D.eventTimes))';
                for r = 1:num_choices %split by choice for plotting PSTHs later
                    popRate(r,t) = mean(D.popRates(D.response==r));
                    popRate_err(r,t) = std(D.popRates(D.response==r));
                end
                
                %extract individual neural activities
                for u=1:numUnits
                    unitID = Units(u);
                    spikeTimes = obj.spikes.shanks(1).units(unitID).spikeTimes;
                    FR_fcn = @(eventTime)(length(find(spikeTimes >= eventTime-0.5*obj.timebinWidth & spikeTimes <= eventTime+0.5*obj.timebinWidth))/obj.timebinWidth);
                    D.firingRate(:,u) = arrayfun(FR_fcn,D.eventTimes + shiftTime);
                end
                
                %Calculate log liks using glmnet
                N=0.4;
                C = cvpartition(D.response,'Kfold',obj.numCVfolds);
                
                p_hats_PR = nan(C.NumObservations,1);
                p_hats_NR = nan(C.NumObservations,1);
                p_hats_baseline = nan(C.NumObservations,1);
                for f = 1:C.NumTestSets
                    trainIdx = C.training(f);
                    testIdx = C.test(f);
                    testID = find(testIdx);
                    
                    Ytrain = D.response(trainIdx);
                    Ytest = D.response(testIdx);
                    
                    %Fit GLM without any neural activity
                    Xtrain = [D.contrast_cond(trainIdx,:).^N];
                    Xtest = [D.contrast_cond(testIdx,:).^N];
                    switch(num_choices)
                        case 2 %AFC
                            fit= glmnet(Xtrain,Ytrain,'binomial');
                            ph = glmnetPredict(fit,Xtest,0.01,'response');
                            ph = [1-ph ph];
                        case 3 %3 choice
                            fit= glmnet(Xtrain,Ytrain,'multinomial');
                            ph = glmnetPredict(fit,Xtest,0.01,'response');
                    end
                    
                    for i=1:length(Ytest)
                        p_hats_baseline(testID(i)) = ph(i,Ytest(i));
                    end
                    
                    %Fit GLM with population rate
                    Xtrain = [D.contrast_cond(trainIdx,:).^N D.popRates(trainIdx)];
                    Xtest = [D.contrast_cond(testIdx,:).^N D.popRates(testIdx)];
                    switch(num_choices)
                        case 2 %AFC
                            fit= glmnet(Xtrain,Ytrain,'binomial');
                            ph = glmnetPredict(fit,Xtest,0.01,'response');
                            ph = [1-ph ph];
                        case 3 %3 choice
                            fit= glmnet(Xtrain,Ytrain,'multinomial');
                            ph = glmnetPredict(fit,Xtest,0.01,'response');
                    end
                    
                    for i=1:length(Ytest)
                        p_hats_PR(testID(i)) = ph(i,Ytest(i));
                    end
                    
                    %Fit GLM with individual neural firing rates
                    Xtrain = [D.contrast_cond(trainIdx,:).^N D.firingRate(trainIdx,:)];
                    Xtest = [D.contrast_cond(testIdx,:).^N D.firingRate(testIdx,:)];
                    
                    switch(num_choices)
                        case 2 %AFC
                            opts.alpha=0; %L2 regularisation
                            opts.penalty_factor = [0 0 ones(1,numUnits)];
                            fit= glmnet(Xtrain,Ytrain,'binomial',glmnetSet(opts));
                            ph = glmnetPredict(fit,Xtest,0.01,'response'); %Predict at lambda=0.01
                            ph = [1-ph ph];
                        case 3 %3 choice
                            opts.alpha=0; %L2 regularisation
                            opts.penalty_factor = [0 0 ones(1,numUnits)];
                            fit= glmnet(Xtrain,Ytrain,'multinomial',glmnetSet(opts));
                            ph = glmnetPredict(fit,Xtest,0.01,'response'); %Predict at lambda=0.01
                    end
                    
                    for i=1:length(Ytest)
                        p_hats_NR(testID(i)) = ph(i,Ytest(i));
                    end
                    
                end
                LL_PR(t) = mean(-log2(p_hats_PR));
                LL_NR(t) = mean(-log2(p_hats_NR));
                LL_baseline(t) = mean(-log2(p_hats_baseline));
                
                
                %Plot the Popmodel loglik
                ax=subplot(2,1,1);
                plot(ax,timeseries,LL_baseline,'kd-',timeseries,LL_PR,'mo:',timeseries,LL_NR,'gs:');
                legend({'Baseline','+PopulationRate','+IndivNeur(L2)'});
                ylabel('Model Likelihood [Bits per trial]');
                xlim(time);
                
                %plot the popRate PSTH split by choice
                ay=subplot(2,1,2);
                cla(ay);
                %                 axes(ay);
                cols = {'b','r','y'};
                hold on;
                for r=1:num_choices
                    errH = popRate(r,:)+popRate_err(r,:);
                    errL = popRate(r,:)-popRate_err(r,:);
                    plotWithErrUL(ay,timeseries,popRate(r,:),[errH; errL],cols{r});
                end
                
                %Plot time markers
                psth_ylim = get(gca,'ylim');
                line([0 0],psth_ylim);
                
                hold off;
                ylabel('Population rate');
                xlim(time);
                %                 ylim(psth_ylim);
                
                drawnow;
            end
            
            if isempty(obj.split.Groups)
                title([epoch ' : all isolated cells']);
                filename = [obj.name '_' epoch];
            else
                title([epoch ' : Cells in division ' num2str(splitFlag)]);
                filename = [obj.name '_' epoch '_cellGroup_' num2str(splitFlag)];
            end
            set(get(gcf,'children'),'fontsize',15);
            figDir = 'B:\figures\GLM+NeuralActivity';
            savefig(fullfile(figDir,[filename '.fig']));
            print(fullfile(figDir,[filename '.pdf' ]),'-dpdf','-painters');
            
        end
    end
    
    methods (Access=private)
        function obj = NoGoResampling(obj)
            disp('Correcting for NoGo response timing...');
            GO_idx = obj.cwLabels.responseMade<3;
            NG_idx = obj.cwLabels.responseMade==3;
            numNG = sum(NG_idx);
            
            GO_delays = obj.cwEvents.responseMade(GO_idx) - obj.cwEvents.goCue(GO_idx);
            obj.cwEvents.responseMade(NG_idx) = obj.cwEvents.goCue(NG_idx) + datasample(GO_delays,numNG);
            
            GO_delays = obj.cwEvents.responseMoveStartTimes(GO_idx) - obj.cwEvents.goCue(GO_idx);
            obj.cwEvents.responseMoveStartTimes(NG_idx) = obj.cwEvents.goCue(NG_idx) + datasample(GO_delays,numNG);
            
            GO_delays =  obj.cwEvents.feedbackStarted(GO_idx) - obj.cwEvents.goCue(GO_idx);
            obj.cwEvents.feedbackStarted(NG_idx) = obj.cwEvents.goCue(NG_idx) + datasample(GO_delays,numNG);
        end
    end
    
end