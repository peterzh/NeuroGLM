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
        ContrastDimensions;
    end
    
    properties (Access=public)
        guess_bpt;
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
                D = struct('contrast_cond',[],'response',[],'repeatNum',[],'feedbackType',[],'RT',[]);
                
                try
                    block = dat.loadBlock(obj.expRef);
                    trials = block.trial;
                    
                    for t=1:block.numCompletedTrials
                        D.contrast_cond(t,:) = trials(t).condition.visCueContrast';
                        D.response(t,1) = trials(t).responseMadeID';
                        D.repeatNum(t,1) = trials(t).condition.repeatNum;
                        D.feedbackType(t,1) = trials(t).feedbackType;
                        D.RT(t,1) = trials(t).responseMadeTime-trials(t).interactiveStartedTime;
                        
%                         if D.RT(t,1) > 4
%                             keyboard;
%                         end
                    end
                catch
                    warning('session data empty');
                end
                
                obj.data = D;
            else
                error('GLM:constructorFail', 'Must pass either an expRef or data struct to GLM constructor');
            end
            
            if ~isempty(obj.data.response)
                
                if any(min(obj.data.contrast_cond,[],2)>0)
                    obj.ContrastDimensions = 2;
                else
                    obj.ContrastDimensions = 1;
                end
                
                tab = tabulate(obj.data.response);
                tab = tab(:,3)/100;
                obj.guess_bpt=sum(tab.*log2(tab));
            else
                obj.ContrastDimensions = NaN;
                obj.guess_bpt = NaN;
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
                case 'C-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf;
                        +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1)   );
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2)  );
                case 'C'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','ScaleR_L','Offset_R','ScaleL_R','ScaleR_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf;
                        +inf +inf +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1) + P(3).*in(:,2));
                    obj.ZR = @(P,in)(P(4) + P(5).*in(:,1) + P(6).*in(:,2));
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
                        +inf +inf +inf +inf 3];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1).^P(5));
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2).^P(5));
                case 'C^NL^NR-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N_L','N_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0;
                        +inf +inf +inf +inf 3 3];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1).^P(5));
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2).^P(6));
                case 'C50'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','ScaleR_L','Offset_R','ScaleL_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf +inf +inf 3 0.8];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(7))./(in(:,1).^P(7) + P(8)^P(7)) + P(3).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)));
                    obj.ZR = @(P,in)(P(4) + P(5).*(in(:,1).^P(7))./(in(:,1).^P(7) + P(8)^P(7)) + P(6).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)));
                case 'C50-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf 3 0.8];
                    
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)));
                    obj.ZR = @(P,in)(P(3) + P(4).*(in(:,2).^P(5))./(in(:,2).^P(5) + P(6)^P(5)));
                case 'C50LR-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N_L','C50_L', 'N_R', 'C50_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001 0 0.001;
                        +inf +inf +inf +inf 3 0.8 3 0.8];
                    
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)));
                    obj.ZR = @(P,in)(P(3) + P(4).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)));    
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
                case 'AFC-C^N-subset-hist-success'
                    obj.parameterLabels = {'Offset','ScaleL_L','ScaleR_R','Hist','N'};
                    obj.parameterBounds = [-inf(1,4) 0;
                        +inf(1,4) inf];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2) D.hist]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(5) + P(3).*in(:,2).^P(5) + P(4).*in(:,3) );
                case 'C^N-subset-Qlearning'
                    obj.parameterLabels = {'ScaleL_L','ScaleR_R','Q_L','Q_R','N'};
                    obj.parameterBounds = [-inf(1,2) 1 1 0.1;
                                            inf(1,2) 1 1 1];
                    obj.Zinput = @(D,B)(behav.GLM_QlearningInputs(D,B));
                    obj.ZL = @(P,in)(P(1).*in(:,1).^P(5) + P(3).*in(:,3) ); %+ P(5).*in(:,5));
                    obj.ZR = @(P,in)(P(2).*in(:,2).^P(5) + P(4).*in(:,4) ); %+ P(5).*in(:,5));
                    obj.parameterStart = @()([-1+2*rand(1,4) 0.5]);
                    
                case 'C^N-subset-2AFC'
                    obj.parameterLabels = {'Offset','ScaleL','ScaleR','N'};
                    obj.parameterBounds = [-inf -inf -inf 0.1;
                        +inf +inf +inf 1];
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(4) + P(3).*in(:,2).^P(4)  );
                    obj.parameterStart = [0 0 0 0.1];
                    
                case 'C^N-subset-Qlearning-2AFC'
                    obj.parameterLabels = {'OffsetL-R','ScaleL-R','QL-R'};
                    obj.parameterBounds = [-inf -inf -inf;
                        +inf +inf +inf];
                    obj.Zinput = @(D)(behav.GLM_QlearningInputs(D));
                    N = 0.4;
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1) + P(3).*in(:,2) );
                    obj.parameterStart = @()([10*rand(1,3)]);
                    
                case 'C^N-subset-Qlearning-noBias-2AFC'
                    obj.parameterLabels = {'ScaleL-R','QL-R'};
                    obj.parameterBounds = [-inf -inf;
                        +inf +inf];
                    obj.Zinput = @(D)(behav.GLM_QlearningInputs(D));
                    N = 0.4;
                    obj.ZL = @(P,in)( P(1).*in(:,1) + P(2).*in(:,2) );
                    obj.parameterStart = @()([10*rand(1,2)]);
                    
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
%             obj.data = obj.getrow(obj.data,6:length(obj.data.response));
            
            %Remove trials with repeats
%             obj.data = obj.getrow(obj.data,obj.data.repeatNum==1);
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
            
            responses = obj.data.response;
            
            if isempty(obj.regularise)
                objective = @(b) (obj.calculateLogLik(b, obj.Zinput(obj.data), responses));
            else
                objective = @(b) (obj.calculateLogLik(b, obj.Zinput(obj.data), responses) + obj.regularise(b));
            end
            
            
            
%             if exist('opti','file')==2 %use opti toolbox if loaded
%                 options = optiset('display','final','solver','NOMAD');
%                 Opt = opti('fun',objective,'bounds',obj.parameterBounds(1,:),obj.parameterBounds(2,:),'x0',obj.parameterStart(),'options',options);
%                 [p,~,exitflag,~] = Opt.solve;
%                 obj.parameterFits = p';
%                 if exitflag < 0
%                     obj.parameterFits = nan(1,length(obj.parameterLabels));
%                 end
%             else
                [obj.parameterFits,~,exitflag] = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
                if ~any(exitflag == [1,2])
                    obj.parameterFits = nan(1,length(obj.parameterLabels));
                end
%             end

        end
        
        function [obj,varargout] = fitCV(obj,varargin)
            %Crossvalidated fitting
            
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            %Trim first 5 trials
%             obj.data = obj.getrow(obj.data,6:length(obj.data.response));
            
            %Remove trials with repeats
%             obj.data = obj.getrow(obj.data,obj.data.repeatNum==1);
            
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
                
%                 if exist('opti','file')==2 %use opti toolbox if loaded
%                     options = optiset('solver','NLOPT');
%                     Opt = opti('fun',objective,'bounds',obj.parameterBounds(1,:),obj.parameterBounds(2,:),'x0',obj.parameterStart(),'options',options);
%                     [p,~,exitflag,~] = Opt.solve;
%                     obj.parameterFits(f,:) = p';
%                     if exitflag < 0
%                         obj.parameterFits = nan(1,length(obj.parameterLabels));
%                     end
%                 else
                    [obj.parameterFits(f,:),~,exitflag] = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
                    if ~any(exitflag == [1,2])
                        obj.parameterFits(f,:) = nan(1,length(obj.parameterLabels));
                    end
%                 end

                
                phat = obj.calculatePhat(obj.parameterFits(f,:), testInputs);
                
                for i = 1:length(testResponse)
                    obj.p_hat(testIdx(i),1) = phat(i,testResponse(i));
                end
                
                LL(f)=mean(-log2(obj.p_hat(testIdx)));
            end
            
            varargout = {LL};
        end
        
        function h = plotData(obj)
                        
            switch(obj.ContrastDimensions)
                case 1
                    contrast1D = obj.data.contrast_cond(:,2) - obj.data.contrast_cond(:,1);
                    uniqueC1D = unique(contrast1D);
                    prop=[];
                    prop_ci=[];
                    for c = 1:length(uniqueC1D)
                        D = obj.getrow(obj.data,contrast1D == uniqueC1D(c));
                        p = sum([D.response==1 D.response==2 D.response==3],1)/length(D.response);
                        
                        [p,pci]=binofit(sum([D.response==1 D.response==2 D.response==3],1),length(D.response),0.05);
                        
                        prop_ci(:,:,c) = pci;
                        
                        prop = [prop;p];
                        %                         bse = sqrt((p.*(1-p)/length(D.response)));
                        %                         binom_se = [binom_se;bse];
                    end
                    
                    
                    
%                     plot(uniqueC1D,prop,'.','MarkerSize',20);
                    
                    err = [prop - squeeze(prop_ci(:,1,:))', squeeze(prop_ci(:,2,:))' - prop];
                    
                    if max(obj.data.response) == 2 %for 2AFC tasks
                        r = 2;
                    else
                        r = 3;
                    end
                    
                    errorbar(repmat(uniqueC1D,1,r),prop(:,1:r),err(:,[1:r]),err(:,[4:(3+r)]),'.','MarkerSize',20)
                    xlabel('contrast');
                    ylabel('P( choice | contrast)');
                    
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
                    
                    titles = {'p( Left | c)','p( Right | c)','p( NoGo | c)'};
                    for i=1:3
                        h(i)=subplot(2,3,i);
                        imagesc(uniqueCR,uniqueCL,prop(:,:,i),[0 1]);
                        set(gca,'YDir','normal','box','off');
                        
%                         xlabel('Contrast right');
%                         ylabel('Contrast left');
                        title(titles{i});
                        axis square;
                        set(gca,'XTick','','YTick',0:0.1:0.5);
                        if i > 1
                            set(gca,'XTick','','ytick','');
                        end
                        
                        if i == 1
                            %                                 xlabel('Contrast right');
                            ylabel('Contrast left');
                        end
                    end
            end
            
            set(gcf,'color','w');
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
%                         titles = {'pred P( left | contrast)','pred P( right | contrast)','pred P( nogo | contrast)'};
                        for i=1:3
                            subplot(2,3,i+3);
                            imagesc(evalCR,evalCL,prop(:,:,i),[0 1]);
                            set(gca,'YDir','normal','box','off','YTick',0:0.1:0.5);
                            
                            
                            if i > 1
                                set(gca,'XTick','','ytick','');
                            end
                            
                            if i == 1
                                xlabel('Contrast right');
                                ylabel('Contrast left');
                            end
%                             xlabel('Contrast right');
%                             ylabel('Contrast left');
%                             title(titles{i});
                            axis square;
                        end
                        
                end
            else
                error('Model not fitted (non-crossvalidated) yet');
            end
        end
        
        function plotPedestal(obj)
            if isempty(obj.parameterFits)
                error('Need to fit model first');
            end
            
            cVals = unique(obj.data.contrast_cond(:));
            prop=nan(length(cVals),length(cVals),3);
            
            for cl = 1:length(cVals)
                for cr = 1:length(cVals)
                    E = obj.getrow(obj.data,obj.data.contrast_cond(:,1) == cVals(cl) & obj.data.contrast_cond(:,2) == cVals(cr));
                    for i=1:3
                        prop(cl,cr,i) = sum(E.response==i)/length(E.response);
                    end
                    pd.propChooseLeft(cl,cr) = prop(cl,cr,1);
                    pd.propChooseRight(cl,cr) = prop(cl,cr,2);
                    pd.propChooseNogo(cl,cr) = prop(cl,cr,3);
                end
            end
            
            cTestVals = min(cVals):0.01:2;
            selectContrast{1} = [reshape(reshape(repmat(cVals, 1, length(cTestVals)), length(cVals), length(cTestVals))', 1, length(cVals)*length(cTestVals)); ...
                repmat(cTestVals, 1, length(cVals))];
            selectContrast{2} = selectContrast{1}([2 1], :);
            
            predictionsSelect{1} = obj.calculatePhat(obj.parameterFits, selectContrast{1}')';
            predictionsSelect{2} = obj.calculatePhat(obj.parameterFits, selectContrast{2}')';
            f3 = figure; %set(f3, 'Position', [  -1896         507        1058         405]);
            
            colors = [        0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250];
            
            for ped = 1:length(cVals)-1
                
                subplot(1, length(cVals)-1, ped)
                
                for r = 1:3
                    plot(-(cVals(ped:end)-cVals(ped)), prop(ped:end,ped,r), 'o','Color',colors(r,:));

                    hold on;
                    plot((cVals(ped:end)-cVals(ped)), prop(ped,ped:end,r), 'o','Color',colors(r,:));

                    
                    
                    plot(-(cTestVals(cTestVals>=cVals(ped))-cVals(ped)), predictionsSelect{2}(r,selectContrast{1}(1,:)==cVals(ped)&selectContrast{1}(2,:)>=cVals(ped)),'Color',colors(r,:));

                    plot((cTestVals(cTestVals>=cVals(ped))-cVals(ped)), predictionsSelect{1}(r,selectContrast{2}(2,:)==cVals(ped)&selectContrast{2}(1,:)>=cVals(ped)),'Color',colors(r,:));

                end
                
                title(['pedestal = ' num2str(cVals(ped))]);
                xlabel('delta C');
                ylim([0 1]);
                xlim([-1 1]);
                set(gca,'box','off');
%                 makepretty
            end
            set(gcf,'color','w');
        end
        
        function plotPredVsActual(obj, varargin)
            
            plotParams.Color = 'r';
            plotParams.ax = gca;
            if ~isempty(varargin)
                plotParams = mergeStructs(varargin{1},plotParams);
            end
            
            switch(obj.ContrastDimensions)
                case 1
                    contrast1D = diff(obj.data.contrast_cond, [], 2);
                    uniqueC1D = unique(contrast1D);
                    nC = length(uniqueC1D);
                    prop=zeros(nC,3);
                    prop_ci=zeros(nC,3,2);
                    for c = 1:length(uniqueC1D)
                        D = obj.getrow(obj.data,contrast1D == uniqueC1D(c));
                        respSum = sum([D.response==1 D.response==2 D.response==3],1);
                        p = respSum/length(D.response);
                        
                        [p,pci]=binofit(respSum,length(D.response),0.05);
                        
                        prop_ci(c,:,:) = pci;
                        
                        prop(c,:) = p;
                    end                                                            
                    
                    if max(obj.data.response) == 2 %for 2AFC tasks
                        rMax = 2;
                    else
                        rMax = 3;
                    end                    

                    evalCon = unique(obj.data.contrast_cond,'rows');
                    evalC1d = evalCon(:,2) - evalCon(:,1);
                    [~,sortIdx]=sort(evalC1d);
                    evalCon = evalCon(sortIdx,:);
                    phat = obj.calculatePhat(obj.parameterFits,evalCon);
                    
                    rSymbols = {'o', '.', 'x'};
                    
                    for c = 1:length(uniqueC1D)
                        for r = 1:rMax
                            plot(plotParams.ax, phat(c,r), prop(c,r), rSymbols{r}, 'Color', plotParams.Color)
                            hold on;                            
                        end
                        for r = 1:rMax
                            plot(plotParams.ax, phat(c,r)*ones(1,2), squeeze(prop_ci(c,r,:)), 'Color', plotParams.Color)
                        end
                    end
                    plot(plotParams.ax, [0 1], [0 1], 'k--');

                    ylabel('actual probability');
                    xlabel('predicted probability');
                    legend({'left resp', 'right resp', 'nogo'}, 'Location', 'Best');

                    axis square
                    box off
                    
                case 2
                    fprintf(1, 'plotPredVsActual not yet implemented for 2D task\n')
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
                pNG = 1 - (pL + pR);
                
                phat = [pL pR pNG];
            end
            
%             if any(phat<0)
%                 keyboard;
%             end
        end
        
        function obj = addData(obj,type)
            %Adds extra data depending on what is required
            
            if strcmp(obj.expRef,'none')
                error('not coded for custom struct inputs');
            end
            
            switch(type)
                case 'lick'
                    block = dat.loadBlock(obj.expRef);
                    trials = block.trial;
                    
                    disp('Loading lick data...');
                    try
                        L=load(dat.expFilePath(obj.expRef, 'Timeline', 'm'));
                        tseries = L.Timeline.rawDAQTimestamps;
                        lickseries = L.Timeline.rawDAQData(:,7);
                        
                        for t=1:block.numCompletedTrials
                            start = trials(t).trialStartedTime;
                            finish = trials(t).trialEndedTime;
                            idx = (start < tseries) & (tseries < finish);
                            
                            obj.data.lickenergy(t,1) = sum(lickseries(idx).^2);
                        end
                    catch
                        warning('No lick data found.. Setting to NaNs');
                        obj.data.lickenergy = nan(block.numCompletedTrials,1);
                    end
                    
                case 'pupil'
                    block = dat.loadBlock(obj.expRef);
                    trials = block.trial;
                    
                    disp('Loading eye movie...');
                    
                    try
                        vidFilename = ['\\zserver\Data\EyeCamera\' obj.expRef(14:end) '\' obj.expRef(1:10) '\' obj.expRef(12) '\eye.mj2'];
                        v = VideoReader(vidFilename);
                        
                        for t=1:block.numCompletedTrials
                            %                             disp(t);
                            start = trials(t).trialStartedTime;
                            finish = trials(t).trialEndedTime;
                            
                            v.CurrentTime = start; a=0; pixCounts=0;
                            while v.CurrentTime < finish
                                frame = 255 - (v.readFrame*5.5);
                                frame(frame>0)=1;
                                frame = frame(50:110,80:180);
                                pixCounts = pixCounts + sqrt(sum(frame(:)));
                                a=a+1;
                            end
                            
                            obj.data.pupilenergy(t,1) = pixCounts/a;
                        end
                    catch
                        warning('No eye data found.. Setting to NaNs');
                        obj.data.pupilenergy = nan(block.numCompletedTrials,1);
                        
                    end
            end
            
            disp('Done!');
            
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