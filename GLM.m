classdef GLM
    properties (Access=public)
        expRef;
        modelString;
        parameterLabels;
        parameterFits;
        parameterBounds;
        parameterStart;
        Zinput;
        ZL;
        ZR;
        regularise;
        data;
        p_hat;
    end
    
    properties (Access=private)
        ContrastDimensions;
    end
    
    methods
        function obj = GLM(inputData)
            if isa(inputData,'struct')
                %If input is a struct containing AT MINIMUM the fields:
                %                 -contrast_cond
                %                 -response
                %                 -repeatNum
                obj.data = inputData;
                obj.expRef = 'none';
                
            elseif isa(inputData,'char')
                %if expRef, then load using the dat package
                obj.expRef = inputData;
                block = dat.loadBlock(obj.expRef);
                trials = block.trial;
                D = struct;
                
                for t=1:block.numCompletedTrials
                    D.contrast_cond(t,:) = trials(t).condition.visCueContrast';
                    D.response(t,1) = trials(t).responseMadeID';
                    D.repeatNum(t,1) = trials(t).condition.repeatNum;
                    D.feedbackType(t,1) = trials(t).feedbackType;
                end

                obj.data = D;
            else
                error('GLM:constructorFail', 'Must pass either an expRef or data struct to GLM constructor');
            end
            
            if any(min(obj.data.contrast_cond,[],2)>0)
                obj.ContrastDimensions = 2;
            else
                obj.ContrastDimensions = 1;
            end
            
        end
        
        function obj = setModel(obj,modelString)
            obj.modelString = modelString;
            obj.parameterFits = [];
            obj.parameterStart = [];
            
            switch(modelString)
                case 'Offset' %Model guesses based on the proportion of responses in the data
                    %used as a baseline to compare other models
                    obj.parameterLabels = {'Offset_L','Offset_R'};
                    obj.parameterBounds = [-inf -inf; +inf +inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1)*ones(length(in(:,1)),1));
                    obj.ZR = @(P,in)(P(2)*ones(length(in(:,1)),1));
                case 'C^N'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','ScaleR_L','Offset_R','ScaleL_R','ScaleR_R','N'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf 0;
                        +inf +inf +inf +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1).^P(7) + P(3).*in(:,2).^P(7));
                    obj.ZR = @(P,in)(P(4) + P(5).*in(:,1).^P(7) + P(6).*in(:,2).^P(7));
                case 'C^N-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0;
                        +inf +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1).^P(5));
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2).^P(5));
                case 'C^NL^NR-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N_L','N_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0;
                        +inf +inf +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1).^P(5));
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2).^P(6));
                case 'C50'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','ScaleR_L','Offset_R','ScaleL_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf 0 0;
                        +inf +inf +inf +inf +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(7))./(in(:,1).^P(7) + P(8)) + P(3).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)));
                    obj.ZR = @(P,in)(P(4) + P(5).*(in(:,1).^P(7))./(in(:,1).^P(7) + P(8)) + P(6).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)));
                case 'C50-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0;
                        +inf +inf +inf +inf +inf +inf];
                    
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)));
                    obj.ZR = @(P,in)(P(3) + P(4).*(in(:,2).^P(5))./(in(:,2).^P(5) + P(6)^P(5)));
                case 'Clogistic-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','k','c0_L','c0_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf -inf;
                        +inf +inf +inf +inf +inf +inf +inf];
                    
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2)./(1 + exp(-P(5).*(in(:,1) - P(6)))) );
                    obj.ZR = @(P,in)(P(3) + P(4)./(1 + exp(-P(5).*(in(:,2) - P(7)))) );
                case 'Supersaturation-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N','C50','Q'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0 0;
                        +inf +inf +inf +inf +inf +inf 2];
                    
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5)*P(7) + P(6)^P(5)*P(7)));
                    obj.ZR = @(P,in)(P(3) + P(4).*(in(:,2).^P(5))./(in(:,2).^P(5)*P(7) + P(6)^P(5)*P(7)));
                case 'AFC'
                    obj.parameterLabels = {'Offset','ScaleL','ScaleR'};
                    obj.parameterBounds = [-inf -inf -inf; +inf +inf +inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2)*in(:,1) + P(3)*in(:,2));
                    obj.ZR = [];
                case 'C^N-subset-hist-response'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','HistL_L','HistR_L','HistNG_L','HistL_R','HistR_R','HistNG_R','N'};
                    obj.parameterBounds = [-inf(1,10) 0;
                                          +inf(1,10) inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2) D.hist(:,1) D.hist(:,2) D.hist(:,3)]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(11) + P(5).*in(:,3) + P(6).*in(:,4) + P(7).*in(:,5) );
                    obj.ZR = @(P,in)( P(3) + P(4).*in(:,2).^P(11) + P(8).*in(:,3) + P(9).*in(:,4) + P(10).*in(:,5) );
                case 'C^N-subset-hist-contrast'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','HistL_L','HistR_L','HistNG_L','HistL_R','HistR_R','HistNG_R','N'};
                    obj.parameterBounds = [-inf(1,10) 0;
                                        +inf(1,10) inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2) D.hist(:,1) D.hist(:,2) D.hist(:,3)]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(11) + P(5).*in(:,3) + P(6).*in(:,4) + P(7).*in(:,5) );
                    obj.ZR = @(P,in)( P(3) + P(4).*in(:,2).^P(11) + P(8).*in(:,3) + P(9).*in(:,4) + P(10).*in(:,5) );
                case 'C^N-subset-hist-success'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','Hist_L','Hist_R','N'};
                    obj.parameterBounds = [-inf(1,6) 0;
                                        +inf(1,6) inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2) D.hist]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(7) + P(5).*in(:,3) );
                    obj.ZR = @(P,in)( P(3) + P(4).*in(:,2).^P(7) + P(6).*in(:,3) );
                otherwise
                    error('Model does not exist');
                    
            end
            
            if isempty(obj.parameterStart)
                obj.parameterStart = zeros(1,length(obj.parameterLabels));
            end

        end
        
        function obj = fit(obj)    
            %Non crossvalidated fitting
            
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
                        
            %Trim first 5 trials
            obj.data = obj.getrow(obj.data,6:length(obj.data.response));
            
            %Remove trials with repeats
            obj.data = obj.getrow(obj.data,obj.data.repeatNum==1);
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
            
            inputs = obj.Zinput(obj.data);
            responses = obj.data.response;
            
            if isempty(obj.regularise)
                objective = @(b) (obj.calculateLogLik(b, inputs, responses));
            else
                objective = @(b) (obj.calculateLogLik(b, inputs, responses) + obj.regularise(b));
            end
            
            if exist('opti','file')==2 %use opti toolbox if loaded
                options = optiset('display','iter','solver','NLOPT');
                Opt = opti('fun',objective,'bounds',obj.parameterBounds(1,:),obj.parameterBounds(2,:),'x0',obj.parameterStart(),'options',options);
                [p,~,exitflag,~] = Opt.solve;
                obj.parameterFits = p';
                if exitflag < 0
                    obj.parameterFits = nan(1,length(obj.parameterLabels));
                end
            else
                [obj.parameterFits,~,exitflag] = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
                if ~any(exitflag == [1,2])
                    obj.parameterFits = nan(1,length(obj.parameterLabels));
                end
            end

        end
        
        function obj = fitCV(obj,varargin)
            %Crossvalidated fitting
            
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            %Trim first 5 trials
            obj.data = obj.getrow(obj.data,6:length(obj.data.response));
            
            %Remove trials with repeats
            obj.data = obj.getrow(obj.data,obj.data.repeatNum==1);
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',2000);
            
            if isempty(varargin)
                C = cvpartition(length(obj.data.response),'LeaveOut');
            else
                C = cvpartition(obj.data.response,'KFold',varargin{1});
            end
            
            obj.parameterFits = nan(C.NumTestSets,length(obj.parameterLabels));
            obj.p_hat = [];
            for f=1:C.NumTestSets
                disp(['Model: ' obj.modelString '. Fold: ' num2str(f) '/' num2str(C.NumTestSets)]);
                trainIdx = find(C.training(f)==1);
                testIdx = find(C.test(f)==1);
                
                inputs = obj.Zinput(obj.data);
                trainInputs = inputs(trainIdx,:);
                testInputs = inputs(testIdx,:);
                
                trainResponses = obj.data.response(trainIdx);
                testResponse = obj.data.response(testIdx);
                
                if isempty(obj.regularise)
                    objective = @(b) ( obj.calculateLogLik(b, trainInputs, trainResponses) );
                else
                    objective = @(b) ( obj.calculateLogLik(b, trainInputs, trainResponses) + obj.regularise(b));
                end
                
                if exist('opti','file')==2 %use opti toolbox if loaded
                    options = optiset('solver','NLOPT');
                    Opt = opti('fun',objective,'bounds',obj.parameterBounds(1,:),obj.parameterBounds(2,:),'x0',obj.parameterStart(),'options',options);
                    [p,~,exitflag,~] = Opt.solve;
                    obj.parameterFits(f,:) = p';
                    if exitflag < 0
                        obj.parameterFits = nan(1,length(obj.parameterLabels));
                    end
                else
                    [obj.parameterFits(f,:),~,exitflag] = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
                    if ~any(exitflag == [1,2])
                        obj.parameterFits(f,:) = nan(1,length(obj.parameterLabels));
                    end
                end

                
                phat = obj.calculatePhat(obj.parameterFits(f,:), testInputs);
                
                for i = 1:length(testResponse)
                    obj.p_hat(testIdx(i),1) = phat(i,testResponse(i));
                end
            end
        end
        
        function h = plotData(obj)
                        
            switch(obj.ContrastDimensions)
                case 1
                    contrast1D = obj.data.contrast_cond(:,2) - obj.data.contrast_cond(:,1);
                    uniqueC1D = unique(contrast1D);
                    prop=[];
                    binom_se=[];
                    for c = 1:length(uniqueC1D)
                        D = obj.getrow(obj.data,contrast1D == uniqueC1D(c));
                        p = sum([D.response==1 D.response==2 D.response==3])/length(D.response);
                        prop = [prop;p];
                        bse = sqrt((p.*(1-p)/length(D.response)));
                        binom_se = [binom_se;bse];
                    end
                    
%                     plot(uniqueC1D,prop,'.','MarkerSize',20);
                    errorbar(repmat(uniqueC1D,1,3),prop,0.5*binom_se,'.','MarkerSize',20)
                    xlabel('Contrast1D');
                    ylabel('% choice');
                    
                    h=gca;
                    
                case 2
                    uniqueCL = unique(obj.data.contrast_cond(:,1));
                    uniqueCR = unique(obj.data.contrast_cond(:,2));
                    prop=nan(length(uniqueCL),length(uniqueCR),3);
                    
                    for cl = 1:length(uniqueCL)
                        for cr = 1:length(uniqueCR)
                            E = obj.getrow(obj.data,obj.data.contrast_cond(:,1) == uniqueCL(cl) & obj.data.contrast_cond(:,2) == uniqueCR(cr));
                            for i=1:3
                                prop(cl,cr,i) = sum(E.response==i)/length(E.response);
                            end
                        end
                    end
                    
                    titles = {'%L','%R','%NG'};
                    for i=1:3
                        h(i)=subplot(2,3,i);
                        imagesc(uniqueCR,uniqueCL,prop(:,:,i),[0 1]);
                        set(gca,'YDir','normal');
                        
                        xlabel('C Right');
                        ylabel('C Left');
                        title(titles{i});
                        axis square;
                    end
            end
        end
        
        function fig = plotFit(obj)
            if size(obj.parameterFits,1)==1
                h=obj.plotData();
                
                switch (obj.ContrastDimensions)
                    case 1
                        hold on;
                                                
                        if ~(strcmp(obj.modelString,'fullContrasts') || strcmp(obj.modelString,'fullContrasts-subset'))
                            maxC = max(max(obj.data.contrast_cond));
                            evalC = [linspace(maxC,0,100)', zeros(100,1);
                                zeros(100,1), linspace(0,maxC,100)'];
                            evalC1d = evalC(:,2) - evalC(:,1);
                            
                            otherInputs = obj.Zinput(obj.data);
                            otherInputs(:,1:2)=[];
                            
                            if isempty(otherInputs)
                                inputs = evalC;
                            else
                                inputs = [evalC, zeros(length(evalC),size(otherInputs,2))];
                            end
                            
                            phat = obj.calculatePhat(obj.parameterFits,inputs);
                            set(gca, 'ColorOrderIndex', 1);
                            plot(h, evalC1d,phat);
                            title(obj.modelString);
                            hold off;
                            h=gca;
                        else
                            evalC = unique(obj.data.contrast_cond,'rows');
                            evalC1d = evalC(:,2) - evalC(:,1);
                            [~,sortIdx]=sort(evalC1d);
                            phat = obj.calculatePhat(obj.parameterFits,evalC);
                            set(gca, 'ColorOrderIndex', 1);
                            plot(h, evalC1d(sortIdx),phat(sortIdx,:),':');
                        end
                        title(obj.modelString);
                        hold off;
                        h=gca;
                        
                    case 2
                        h=obj.plotData;
                        fig=get(h(1),'Parent');
                        
                        evalCL = linspace(0,max(obj.data.contrast_cond(:,1)),100);
                        evalCR = linspace(0,max(obj.data.contrast_cond(:,1)),100);
                        prop=nan(length(evalCL),length(evalCR),3);
                        
                        for cl = 1:length(evalCL)
                            for cr = 1:length(evalCR)
                                p = obj.calculatePhat(obj.parameterFits,[evalCL(cl) evalCR(cr)]);
                                for i=1:3
                                    prop(cl,cr,i) = p(i);
                                end
                            end
                        end
                        
                        figure(fig);
                        titles = {'Pred %L','Pred %R','Pred %NG'};
                        for i=1:3
                            subplot(2,3,i+3);
                            imagesc(evalCR,evalCL,prop(:,:,i),[0 1]);
                            set(gca,'YDir','normal');
                            xlabel('C Right');
                            ylabel('C Left');
                            title(titles{i});
                            axis square;
                        end
                        
                end
            else
                error('Model not fitted (non-crossvalidated) yet');
            end
        end
        
        function h = plotParams(obj)
            if size(obj.parameterFits,1)==1
                bar(obj.parameterFits);
                set(gca,'XTickLabel',obj.parameterLabels,'XTick',1:numel(obj.parameterLabels));
                title(obj.modelString);
                h=gca;
            end
        end
        
        function phat = calculatePhat(obj,testParams,inputs)
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end

            if isempty(obj.ZR) %if a AFC task then no ZR is defined, only pL vs pR
                zl = obj.ZL(testParams,inputs);
                pL = exp(zl)./(1+exp(zl));
                pR = 1 - pL;
                N = length(pL);
                phat = [pL pR zeros(N,1)];
            else 
                zl = obj.ZL(testParams,inputs);
                zr = obj.ZR(testParams,inputs);
                pL = exp(zl)./(1+exp(zl)+exp(zr));
                pR = exp(zr)./(1+exp(zl)+exp(zr));
                pNG = 1 - pL - pR;
                
                phat = [pL pR pNG];
            end
        end
    end
    
    methods (Access= {?GLM})        
        function logLik = calculateLogLik(obj,testParams, inputs, responses)
            phat = obj.calculatePhat(testParams, inputs);
            logLik = -sum(log( phat(sub2ind(size(phat), [1:length(responses)]', responses)) ));
        end
        
        function row = getrow(~,D,numrow)
            % Version 1.0 9/18/03
            % by Joern Diedrichsen
            % http://www.icn.ucl.ac.uk/motorcontrol/toolboxes/toolbox_util.htm
            
            if (~isstruct(D))
                error('D must be a struct');
            end;
            
            field = fieldnames(D);
            row=[];
            for f=1:length(field)
                F = getfield(D,field{f});
                if iscell(F)
                    row = setfield(row,field{f},F(numrow,:));
                else
                    row = setfield(row,field{f},F(numrow,:));
                end
            end
            
        end
    end
end