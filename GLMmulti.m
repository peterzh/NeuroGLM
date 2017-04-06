classdef GLMmulti
    properties
        expt;
        names;
        expRefs;
        data;
        model;
        fitData;
        fitLambda = 0.0001;
        guess_bpt = [];
    end
    
    properties (Access=public)
        significance_flag = 0;
        moreData = 0;
        nestedModel_flag = 0;
    end
    
    methods
        function obj = GLMmulti(expt,nameList)
            %default names
            if isempty(nameList)
                nameList = {'Spemann','Murphy','Morgan','Whipple',...
                            'Chomsky','Laveran','Eijkman','Erlanger',...
                            'Gasser','Fleming','Chain','Muller',...
                            'Radnitz','Cori','Moniz','M140528_NS1',...
                            'M140617_NS1','M140701_NS1'};
            end
            
            disp('Collating expRefs...');
            obj.names = nameList;
            obj.expt = expt;
            numSubjects = length(obj.names);
            obj.expRefs = cell(1,numSubjects);
            for n = 1:numSubjects
                obj.expRefs{n} = obj.getExpRefs(obj.names{n},expt);
            end
            
            empty = cellfun(@isempty,obj.expRefs);
            if any(empty)
                error(['No data found for ' obj.names{empty} ' please exclude']);
            end
            
            disp('Loading data...');
            obj.data = cell(1,numSubjects);
            for n = 1:numSubjects
                D = struct;
                for session = 1:length(obj.expRefs{n})
                    eRef = obj.expRefs{n}{session};
                    d = struct('stimulus',[],'response',[],'repeatNum',[],'feedbackType',[],'RT',[]);
                    
                    %if its a laser session, then include only the
                    %nonlaser trials
                    if exist(dat.expFilePath(eRef,'lasermanip','m'), 'file') == 2
                        l = laserGLM(eRef);
                        d = l.data;
                        d.stimulus = d.contrast_cond;
                        d = getrow(d,d.laserIdx==0);
                        d = rmfield(d,{'contrast_cond','laserIdx','laser'});
                    else
                        block = dat.loadBlock(eRef);
                        trials = block.trial;
                        
                        hist=nan(block.numCompletedTrials,12);
                        for t=1:block.numCompletedTrials
                            d.stimulus(t,:) = trials(t).condition.visCueContrast';
                            d.response(t,1) = trials(t).responseMadeID';
                            d.repeatNum(t,1) = trials(t).condition.repeatNum;
                            d.feedbackType(t,1) = trials(t).feedbackType;
                            d.RT(t,1) = trials(t).responseMadeTime-trials(t).interactiveStartedTime;
                            
                            %                         if t > 5
                            %                             hist(t,:) = [NaN d.response(t-5:t-1)' NaN d.feedbackType(t-5:t-1)'];
                            %                         end
                        end
                    end
                    
%                     d.stimulus = [d.stimulus hist];
                    
                    d = structfun(@(x)(x(6:(end-14),:)),d,'uni',0); %trim first 5 trials and last 15
                    d.sessionID = ones(length(d.response),1)*session;
                    D = addstruct(D,d);
                end
                
                D = getrow(D,D.repeatNum==1);
                obj.data{n} = D;
                
                tab = tabulate(D.response);
                tab = tab(:,3)/100;
                obj.guess_bpt(n)=sum(tab.*log2(tab));
            end
            
            %Create new subject which is the combined data of all subjects
            obj.names = [obj.names '<<ALLSUBJ>>'];
            obj.expRefs = [obj.expRefs {cat(1,obj.expRefs{:})} ];
            obj.data{numSubjects+1} = obj.data{1};
            prevMaxSID = 0;
            for n = 2:numSubjects
                prevMaxSID = prevMaxSID + max(obj.data{n-1}.sessionID);
                addDat = obj.data{n};
                addDat.sessionID = addDat.sessionID+prevMaxSID;
                obj.data{numSubjects+1} = addstruct(obj.data{numSubjects+1},addDat);
            end
            tab = tabulate(obj.data{numSubjects + 1}.response);
            tab = tab(:,3)/100;
            obj.guess_bpt(numSubjects + 1)=sum(tab.*log2(tab));
            
        end
        
        function m = getModel2(~,model)
            %More efficient model definitions
            m.name = model;
            switch(model)
                case 'BO_C^N'
                    m.cfn = @(c,n)(c.^n);
                    m.ZL = @(p,c)( p(1) + 10*p(3).*m.cfn(c(:,1),p(5)) );
                    m.ZR = @(p,c)( p(2) + 10*p(4).*m.cfn(c(:,2),p(5)) );
                    m.LB = [-inf -inf 0 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 1 0]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','sL','sR','nL','nR'};
                    
                case 'BC_C^N'
                    m.cfn = @(c,n)(c.^n);
                    m.ZL = @(p,c)( 10*p(3).*(m.cfn(c(:,1),p(5)) + p(1)) );
                    m.ZR = @(p,c)( 10*p(4).*(m.cfn(c(:,2),p(5)) + p(2)) );
                    m.LB = [-inf -inf 0 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 1 0]; %only used when doing non-laser fit
                    m.pLabels = {'bL','bR','sL','sR','nL','nR'};
                    
                case 'BO_NC50'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p,c)( p(1) + 10*p(3).*m.cfn(c(:,1),p(5),p(7)) );
                    m.ZR = @(p,c)( p(2) + 10*p(4).*m.cfn(c(:,2),p(5),p(7)) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','sL','sR','nL','nR','c50L','c50R'};
                    
                case 'BO_NC50_DISTORTION'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p,c)( p(1) + 10*p(3).*( m.cfn(c(:,1),p(5),p(7)) - m.cfn(c(:,2),p(5),p(7)).*(m.cfn(c(:,1),p(5),p(7)) - p(9))*1/3 ) );
                    m.ZR = @(p,c)( p(2) + 10*p(4).*( m.cfn(c(:,2),p(5),p(7)) - m.cfn(c(:,1),p(5),p(7)).*(m.cfn(c(:,2),p(5),p(7)) - p(10))*1/3 ) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0 1 1]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','sL','sR','nL','nR','c50L','c50R','phiL','phiR'};
                    
                case 'BO_NC50_NORMALISATION'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p,c)( p(1) + 10*p(9).*m.cfn(c(:,1),p(5),p(7)) + 10*p(3).*m.cfn(c(:,1),p(5),p(7))./(m.cfn(c(:,1),p(5),p(7)) + m.cfn(c(:,2),p(5),p(7))) );
                    m.ZR = @(p,c)( p(2) + 10*p(10).*m.cfn(c(:,2),p(5),p(7)) + 10*p(4).*m.cfn(c(:,2),p(5),p(7))./(m.cfn(c(:,1),p(5),p(7)) + m.cfn(c(:,2),p(5),p(7))) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0 +inf +inf]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','snL','snR','nL','nR','c50L','c50R','sL','sR'};
                    
                    
                case 'BO_NC50_2SIDE'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p,c)( p(1) + 10*p(3).*m.cfn(c(:,1),p(5),p(7)) + 10*p(5).*m.cfn(c(:,2),p(5),p(7)) );
                    m.ZR = @(p,c)( p(2) + 10*p(4).*m.cfn(c(:,2),p(5),p(7)) + 10*p(6).*m.cfn(c(:,1),p(5),p(7)) );
                    m.LB = [-inf -inf 0 0 0 0 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','sL','sR','sL.','sR.','nL','nR','c50L','c50R'};
                    
                    
                case 'BC_NC50'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    %                     m.cfn = @(c,varargin)((c.^varargin{1})./(c.^varargin{1} + varargin{2}.^varargin{1}));
                    m.ZL = @(p,c)( 10*p(3).*(m.cfn(c(:,1),p(5),p(7)) + p(1)) );
                    m.ZR = @(p,c)( 10*p(4).*(m.cfn(c(:,2),p(5),p(7)) + p(2)) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'bL','bR','sL','sR','nL','nR','c50L','c50R'};
                    
                    
                case 'BO_NC50_NESTED'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p,c)( p(1) + p(3).*max([m.cfn(c(:,1),p(5),p(7)) m.cfn(c(:,2),p(5),p(7))],[],2) ); %actually ZGnG
                    m.ZR = @(p,c)( p(2) + p(4).*diff([m.cfn(c(:,1),p(5),p(7)) m.cfn(c(:,2),p(5),p(7))],[],2) ); %actually ZLR
                    m.LB = [-inf -inf -inf -inf 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'oGO','oLR','sGO','sLR','nL','nR','c50L','c50R'};
                    
                otherwise
                    error(['Model ' model ' does not exist']);
            end
            
        end
        
        function obj = fit(obj,model)
            %             f=figure;
            obj.model = obj.getModel2(model);
            obj.fitData.perSession = 1;
            
            if strcmp(model(end-5:end),'NESTED')
                obj.nestedModel_flag = 1;
            else
                obj.nestedModel_flag = 0;
            end
            
            numSubjects = length(obj.names);
            obj.fitData.params = cell(1,numSubjects);
            
            for n = 1:numSubjects
                
                sessions = unique(obj.data{n}.sessionID);
                if obj.fitData.perSession==1
                    params = nan(length(sessions),length(obj.model.pLabels));
                    for s = 1:length(sessions)
                        cont = obj.data{n}.stimulus(obj.data{n}.sessionID==sessions(s),1:2);
                        resp = obj.data{n}.response(obj.data{n}.sessionID==sessions(s));
                        
                        params(s,:) = obj.util_optim(obj.model.ZL,obj.model.ZR,cont,resp,obj.model.LB,obj.model.UB);
                    end
                else
                    params = nan(1,length(obj.model.pLabels));
                    cont = obj.data{n}.stimulus;
                    resp = obj.data{n}.response;
                    
                    params(1,:) = obj.util_optim(obj.model.ZL,obj.model.ZR,cont,resp,obj.model.LB,obj.model.UB);
                    params = repmat(params,length(sessions),1);
                end
                
                
                %If the base model does not fit shape params on each
                %side, make sure to copy the parameter used for each side
                params(:,obj.model.LB==0 & obj.model.UB==0) = params(:,find(obj.model.LB==0 & obj.model.UB==0) - 1);
                
                obj.fitData.params{n} = params;
                
            end
        end
        
        
        function obj = bootstrap(obj) %OBSOLETE
            %Fits the non-laser model for one session with resampling, to
            %see parameter tradeoff
            %Just use GLMNET fitting as that's quicker to code
            numSubjects = length(obj.names);
            
            figure;
            for n = 1:numSubjects
                [ZL_nL,ZR_nL,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                nL = obj.data{n}.laserIdx==0; %& obj.data{n}.sessionID==1;
                tab = tabulate(obj.data{n}.sessionID);
                numTrialsToSample = round(mean(tab(:,2)));
                p_nL = [];
                
                for iter = 1:200
                    %                     nL_bootstrap = randsample(find(nL),sum(nL),true);
                    nL_bootstrap = randsample(find(nL),numTrialsToSample,true);
                    cont = obj.data{n}.stimulus(nL_bootstrap,:);
                    resp = obj.data{n}.response(nL_bootstrap);
                    
                    offset = zeros(length(resp),4);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.laserLambda*sum(PARAMETERS.^2) );
                    
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                            p_nL(iter,:) = Opt.solve';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                            p_nL(iter,:) = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                    end
                end
                
                if numP == 2
                    subplot(1,numSubjects,n); %bL vs bR
                    plot(p_nL(:,1),p_nL(:,2),'k.'); hold on;
                    plot(obj.fitData.nonLaserParams{n}(:,1),obj.fitData.nonLaserParams{n}(:,2),'ro');
                    xlabel('B_L'); ylabel('B_R'); title(obj.names{n});
                elseif numP == 4
                    subplot(2,numSubjects,n); %bL vs sL
                    plot(p_nL(:,1),p_nL(:,3),'k.'); hold on;
                    plot(obj.fitData.nonLaserParams{n}(:,1),obj.fitData.nonLaserParams{n}(:,3),'ro');
                    xlabel('B_L'); ylabel('S_L'); title(obj.names{n});
                    
                    subplot(2,numSubjects,numSubjects+n); %bL vs sL
                    plot(p_nL(:,2),p_nL(:,4),'k.'); hold on;
                    plot(obj.fitData.nonLaserParams{n}(:,2),obj.fitData.nonLaserParams{n}(:,4),'ro');
                    xlabel('B_R'); ylabel('S_R');
                end
                drawnow;
            end
            
        end
        
        function obj = removeSubj(obj)
            %Deletes all individual subject's data, keeping only pooled data
            numSubjects = length(obj.names);
            obj.expRefs(1:numSubjects-1)=[];
            obj.data(1:numSubjects-1)=[];
            obj.guess_bpt(1:numSubjects-1)=[];
            obj.names(1:numSubjects-1)=[];
        end
        
        function crossval(obj)
            %For each mouse, assess different nonlaser models to account
            %for behavioural data (only on nonlaser trials).
            
            
            numFolds = 20;
            models = {'B','B & C','B & C^N','B & NR'};%,'nes B','nes B & C','nes B & C^N','nes B & NR'};
            numSubjects = length(obj.names);
            figure('name','model comparison','color','w');
            cfn = @(c,n,c50) ( (c.^n)./(c50.^n + c.^n) );
            
            %             axes; hold on;
            for n = 1:numSubjects
                numSessions = max(obj.data{n}.sessionID);
                
                logLik_perSession = nan(numSessions,length(models));
                for session=1:numSessions
                    D = getrow(obj.data{n},obj.data{n}.laserIdx==0 & obj.data{n}.sessionID==session);
                    D.phat = nan(length(D.response),length(models));
                    
                    cv=cvpartition(D.response,'kfold',numFolds);
                    for fold = 1:cv.NumTestSets
                        trainC = D.stimulus(cv.training(fold),1:2);
                        trainR = D.response(cv.training(fold));
                        testC = D.stimulus(cv.test(fold),1:2);
                        testR = D.response(cv.test(fold));
                        
                        for m = 1:length(models)
                            disp(m);
                            switch(m)
                                case 1 %mnr B
                                    ZL = @(p_offset,p,c)( repmat(p(1),size(c,1),1) );
                                    ZR = @(p_offset,p,c)( repmat(p(2),size(c,1),1) );
                                    LB = [-inf -inf ]; %only used when doing non-laser fit
                                    UB = [+inf +inf ]; %only used when doing non-laser fit
                                    pLabels = {'bL','bR'};
                                    obj.nestedModel_flag = 0;
                                case 2 %mnr B & C
                                    ZL = @(p_offset,p,c)( p(1) + 10*p(3).*c(:,1) );
                                    ZR = @(p_offset,p,c)( p(2) + 10*p(4).*c(:,2) );
                                    LB = [-inf -inf 0 0 ]; %only used when doing non-laser fit
                                    UB = [+inf +inf +inf +inf ]; %only used when doing non-laser fit
                                    pLabels = {'bL','bR','sL','sR'};
                                    obj.nestedModel_flag = 0;
                                case 3 %mnr B & C^N
                                    ZL = @(p_offset,p,c)( p(1) + 10*p(3).*c(:,1).^p(5) );
                                    ZR = @(p_offset,p,c)( p(2) + 10*p(4).*c(:,2).^p(5) );
                                    LB = [-inf -inf 0 0 0]; %only used when doing non-laser fit
                                    UB = [+inf +inf +inf +inf 1]; %only used when doing non-laser fit
                                    pLabels = {'bL','bR','sL','sR','n'};
                                    obj.nestedModel_flag = 0;
                                    %                                 case 4 %mnr B & C^NLR
                                    %                                     ZL = @(p_offset,p,c)( p(1) + 10*p(3).*c(:,1).^p(5) );
                                    %                                     ZR = @(p_offset,p,c)( p(2) + 10*p(4).*c(:,2).^p(6) );
                                    %                                     LB = [-inf -inf 0 0 0 0]; %only used when doing non-laser fit
                                    %                                     UB = [+inf +inf +inf +inf 1 1]; %only used when doing non-laser fit
                                    %                                     pLabels = {'bL','bR','sL','sR','nL','nR'};
                                    %                                     obj.nestedModel_flag = 0;
                                case 4 %mnr B & NR
                                    ZL = @(p_offset,p,c)( p(1) + 10*p(3).*cfn(c(:,1),p(5),p(6)) );
                                    ZR = @(p_offset,p,c)( p(2) + 10*p(4).*cfn(c(:,2),p(5),p(6))  );
                                    LB = [-inf -inf 0 0 0.3 0.001]; %only used when doing non-laser fit
                                    UB = [+inf +inf +inf +inf 20  3]; %only used when doing non-laser fit
                                    pLabels = {'bL','bR','sL','sR','n','c50'};
                                    obj.nestedModel_flag = 0;
                                    %                                 case 6 %mnr B & C50LR
                                    %                                     ZL = @(p_offset,p,c)( p(1) + 10*p(3).*cfn(c(:,1),p(5),p(7)) );
                                    %                                     ZR = @(p_offset,p,c)( p(2) + 10*p(4).*cfn(c(:,2),p(6),p(8)) );
                                    %                                     LB = [-inf -inf 0 0 0.3 0.3 0.001 0.001]; %only used when doing non-laser fit
                                    %                                     UB = [+inf +inf +inf +inf 20 20 3 3]; %only used when doing non-laser fit
                                    %                                     pLabels = {'bL','bR','sL','sR','nL','nR','c50L','c50R'};
                                    %                                     obj.nestedModel_flag = 0;
                                    
                                case 5 %nested B
                                    ZL = @(p_offset,p,c)( p(1) ); %actually pGO/pNG
                                    ZR = @(p_offset,p,c)( p(2) ); %actually pL/pR given GO
                                    LB = [-inf -inf ]; %only used when doing non-laser fit
                                    UB = [+inf +inf ]; %only used when doing non-laser fit
                                    pLabels = {'bGO','bLR'};
                                    obj.nestedModel_flag = 1;
                                    
                                case 6 %nested B & C
                                    ZL = @(p_offset,p,c)( p(1) + 10*p(3).*max(c,[],2)  ); %actually pGO/pNG
                                    ZR = @(p_offset,p,c)( p(2) + 10*p(4).*diff(c,[],2) ); %actually pL/pR given GO
                                    LB = [-inf -inf -inf -inf]; %only used when doing non-laser fit
                                    UB = [+inf +inf +inf +inf]; %only used when doing non-laser fit
                                    pLabels = {'bGO','bLR','sGO','sLR'};
                                    obj.nestedModel_flag = 1;
                                case 7 %nested B & C^N
                                    ZL = @(p_offset,p,c)( p(1) + 10*p(3).*(max(c.^p(5),[],2)  ) ); %actually pGO/pNG
                                    ZR = @(p_offset,p,c)( p(2) + 10*p(4).*(diff(c.^p(5),[],2) ) ); %actually pL/pR given GO
                                    LB = [-inf -inf -inf -inf 0]; %only used when doing non-laser fit
                                    UB = [+inf +inf +inf +inf 1]; %only used when doing non-laser fit
                                    pLabels = {'bGO','bLR','sGO','sLR','n'};
                                    obj.nestedModel_flag = 1;
                                case 8 %nested B & NR
                                    ZL = @(p_offset,p,c)( p(1) + 10*p(3).*(max(cfn(c,p(5),p(6)),[],2)  ) ); %actually pGO/pNG
                                    ZR = @(p_offset,p,c)( p(2) + 10*p(4).*(diff(cfn(c,p(5),p(6)),[],2) ) ); %actually pL/pR given GO
                                    LB = [-inf -inf -inf -inf 0.3 0]; %only used when doing non-laser fit
                                    UB = [+inf +inf +inf +inf 20 1]; %only used when doing non-laser fit
                                    pLabels = {'bGO','bLR','sGO','sLR','n','c50'};
                                    obj.nestedModel_flag = 1;
                            end
                            
                            [PARAM,exitflag] = obj.util_optim(ZL,ZR,[],trainC,trainR,LB,UB);
                            if exitflag == 1 %only test if model succeeded in fitting
                                p = obj.calculatePhat(ZL,ZR,[],PARAM,testC);
                                D.phat(cv.test(fold),m) = p(:,1).*(testR==1) + p(:,2).*(testR==2) + p(:,3).*(testR==3);
                            end
                        end
                    end
                    
                    logLik_perSession(session,:) = nanmean(log2(D.phat),1);
                end
                
                subplot(2,round(numSubjects/2),n);
                ax=notBoxPlot(logLik_perSession);
                for i = 1:length(models)
                    ax(i).data.MarkerSize=3;
                    ax(i).data.MarkerFaceColor=[1 1 1]*0.2;
                    ax(i).data.MarkerEdgeColor=[1 1 1]*0.2;
                end
                
                title(obj.names{n});
                
                if n == 1
                    ylabel('CV Log_2 likelihood');
                end
                set(gca,'xtick',1:(length(models)),'XTickLabel',models,'XTickLabelRotation',45,'box','off');
                drawnow;
                
                %                 keyboard;
            end
            
        end
        
        function plotFitOneSession(obj,n,sess)
            if isempty(sess)
                p_nL = median(obj.fitData.params{n},1);
                if obj.fitData.perSession==1
                    warning('Model fits cannot be visualised well here as this requires averaging nonlaser params');
                end
                cont = obj.data{n}.stimulus;
                resp = obj.data{n}.response;
            else
                p_nL = obj.fitData.params{n}(sess,:);
                idx = obj.data{n}.sessionID==sess;
                cont = obj.data{n}.stimulus(idx,:);
                resp = obj.data{n}.response(idx);
            end

            
            cVals = unique(cont(:));
            numPedestals = length(cVals);
            
            figure('name',[obj.names{n} ' session ' num2str(sess)],'color','w');
            %psych curve: 2D representation
            p = nan(length(cVals),length(cVals),3);
            p_pred = nan(size(p));
            for cl = 1:length(cVals)
                for cr = 1:length(cVals)
                    r = resp(cont(:,1) == cVals(cl) & cont(:,2) == cVals(cr));
                    p(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                    p_pred(cl,cr,:) = obj.calculatePhat(obj.model.ZL,obj.model.ZR,p_nL,[cVals(cl) cVals(cr)]);
                end
            end
            
                        labels = {'pL','pR','pNG'};
                        for r = 1:3
                            subplot(3,4,r);
                            imagesc(p(:,:,r)); set(gca,'ydir','normal'); title(['Actual ' labels{r}]);
                            caxis([0 1]);
                            set(gca,'xtick','','ytick',1:length(cVals),'yticklabels',cVals);
                            xlabel('CR'); ylabel('CL');
                            subplot(3,4,r+4);
                            imagesc(p_pred(:,:,r)); set(gca,'ydir','normal'); title(['Pred ' labels{r}]);
                            caxis([0 1]);
                            set(gca,'xtick','','ytick',1:length(cVals),'yticklabels',cVals);
                        end
                        axis(get(gcf,'children'),'square');
%                         set(get(gcf,'children'),'ytick',cVals,'xtick',cVals);
            
            
            %pych curve: pedestal representation
            
            
%             cols = [0 0.4470 0.7410;
%                 0.8500 0.3250 0.0980;
%                 0.4940    0.1840    0.5560];
%             for ped = 1:numPedestals
%                 subplot(numPedestals,1,ped); hold on;
%                 set(gca,'colororder',cols);
%                 %Plot actual datapoints
%                 ped_idx = min(cont,[],2)==cVals(ped);
%                 ped_c_diff = diff(cont(ped_idx,:),[],2);
%                 ped_r = resp(ped_idx);
%                 uC = unique(ped_c_diff);
%                 ph=[];
%                 for c = 1:length(uC)
%                     r = ped_r(ped_c_diff==uC(c));
%                     [ph(c,:),pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
%                     for ch=1:3
%                         l(1)=line([1 1]*uC(c),pci(ch,:));
%                         %                         l(2)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,1));
%                         %                         l(3)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,2));
%                         set(l,'Color',cols(ch,:),'Linewidth',0.5);
%                         %                             l.Color = cols{ch};
%                         %                             l.LineWidth=1;
%                     end
%                 end
%                 set(gca,'ColorOrderIndex',1);
%                 plot(uC,ph,'.','markersize',15);
%                 
%                 %Plot predictions
%                 testCont = [linspace(max(abs(uC))+0.1,0,100)' zeros(100,1); zeros(100,1) linspace(0,max(abs(uC))+0.1,100)'] + cVals(ped);
%                 p_hat = obj.calculatePhat(obj.model.ZL,obj.model.ZR,p_nL,testCont);
%                 set(gca,'ColorOrderIndex',1);
%                 plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
%                 xlim([-1 1]*(max(cVals)+0.1)); ylim([0 1]);
%                 
%                 
%                 title(['pedestal= ' num2str(cVals(ped))]);
%                 axis square;
%                 if ped == numPedestals
%                     ylabel('Choice probability');
%                     xlabel('\Delta contrast from pedestal');
%                 end
%                 if ped < numPedestals
%                     set(gca,'xtick','','xcolor','w');
%                 end
%                 
%             end
%             
            %             %chronometry: RTs
            %             rt = obj.data{n}.RT(idx);
            %             %             keyboard;
            %             CL = (cont(:,1) > cont(:,2)) & resp==1;
            %             CR = (cont(:,1) < cont(:,2)) & resp==2;
            %             hist1=subplot(numPedestals+2,1,numPedestals+1);
            % %             hist(rt(CL));
            %             distributionPlot(rt(CL),'histOpt',1,'histOri','right','xyOri','flipped','showMM',0); axis square;
            %             xlabel('RT correct left [sec]'); set(gca,'ytick','','ycolor','w');
            %             hist2=subplot(numPedestals+2,1,numPedestals+2);
            % %             hist(rt(CR));
            %             distributionPlot(rt(CR),'histOpt',1,'histOri','right','xyOri','flipped','showMM',0); axis square;
            %             xlabel('RT correct right [sec]'); set(gca,'ytick','','ycolor','w');
            %             linkaxes([hist1 hist2],'xy');
            %             set(get(gcf,'children'),'box','off');
            %             xlim([0 1]);
        end
        
        function plotFit(obj)
            numSubjects = length(obj.names);
            
            %Plot session-by-session bias and sensitivity parameters
            figure('name','Session by session non-laser bias and sensitivity changes','color','w');
            i=1;
            for p = [1 2 3 4 5 7]
                h(i)=subplot(3,2,i); hold on;
                
                for n = 1:(length(obj.names)-1)
                    numSessions = length(obj.expRefs{n});
                    boxplot(obj.fitData.params{n}(:,p),'positions',ones(numSessions,1)*n,'plotstyle','compact');
                end
                %                 set(gca,'XTickLabel',obj.names(1:end-1),'xtick',1:(length(obj.names)-1),'xticklabelrotation',45);
                set(gca,'xtick','','xcolor','w');
                set(gca,'box','off');
                
                if i>1
                    %                     set(gca,'ytick','','ycolor','w');
                end
                xlim([0.5 5.5]);
                i=i+1;
                title(obj.model.pLabels{p})
            end
            %             linkaxes(h,'y'); ylim([-6 3]);
            
            
            %Plot predictions against actual psychometric data.
            cols = {[0 0.4470 0.7410],...
                [0.8500 0.3250 0.0980],...
                [0.9290 0.6940 0.1250]};
            
            %             pCorrect = cell(1,numSubjects);
            for n = 1:numSubjects
                figure('name',[obj.names{n} ' non-laser predictions against actual psychometric data']);
                %                 [ZL_nL,ZR_nL,LB] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                
                numSessions = max(obj.data{n}.sessionID);
                for s = 1:numSessions
                    p_nL = obj.fitData.params{n}(s,:);
                    
                    cont = obj.data{n}.stimulus(idx,1:2);
                    resp = obj.data{n}.response(idx);
                    
                    cVals = unique(cont(:));
                    numPedestals = length(cVals)-1;
                    
                    for ped = 1:numPedestals
                        subplot(numPedestals,numSessions,(ped-1)*numSessions + s); hold on;
                        %                         if ped == 1
                        %                             title(obj.names{n})
                        %                         end
                        
                        %Plot actual datapoints
                        ped_idx = min(cont,[],2)==cVals(ped);
                        ped_c_diff = diff(cont(ped_idx,1:2),[],2);
                        ped_r = resp(ped_idx);
                        uC = unique(ped_c_diff);
                        ph=[];
                        for c = 1:length(uC)
                            r = ped_r(ped_c_diff==uC(c));
                            [ph(c,:),pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
                            for ch=1:3
                                l(1)=line([1 1]*uC(c),pci(ch,:));
                                l(2)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,1));
                                l(3)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,2));
                                set(l,'Color',cols{ch},'Linewidth',0.5);
                                %                             l.Color = cols{ch};
                                %                             l.LineWidth=1;
                            end
                        end
                        set(gca,'ColorOrderIndex',1);
                        plot(uC,ph,'.','markersize',15);
                        
                        %Plot predictions
                        testCont = [linspace(max(abs(uC))+0.1,0,100)' zeros(100,1); zeros(100,1) linspace(0,max(abs(uC))+0.1,100)'] + cVals(ped);
                        p_hat = obj.calculatePhat(obj.model.ZL,obj.model.ZR,[],p_nL,testCont);
                        set(gca,'ColorOrderIndex',1);
                        plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
                        xlim([-1 1]*(max(cVals)+0.1)); ylim([0 1]);
                        
                        if s == 1
                            ylabel(['ped: ' num2str(cVals(ped))]);
                        else
                            set(gca,'ytick','','ycolor','w');
                        end
                        
                        if ped ~= numPedestals
                            set(gca,'xtick','','xcolor','w');
                        end
                        
                    end
                    
                    %                     % calculate % correct classification in the dataset
                    %                     % using only easy CL and CR trials
                    %                     easyContrasts = (cont(:,1)==0.54 & cont(:,2)==0) | (cont(:,2)==0.54 & cont(:,1)==0);
                    %                     p_hat = obj.calculatePhat(ZL_nL,ZR_nL,[],p_nL,cont(easyContrasts,:));
                    %                     re = resp(easyContrasts);
                    %                     classify = nan(length(re),1);
                    %                     for t = 1:length(re)
                    %                         [~,rhat] = max(p_hat(t,:));
                    %
                    %                         if rhat == re(t)
                    %                             classify(t)=1;
                    %                         else
                    %                             classify(t)=0;
                    %                         end
                    %                     end
                    %
                    %                     pCorrect{n}(s) = mean(classify);
                    %                     keyboard;
                end
                set(gcf,'color','w');
            end
            
            %             figure('name','model classification accuracy','color','w');
            %             axes; hold on;
            %             for n = 1:(numSubjects-1)
            % %                 subplot(1,numSubjects-1,n);
            %                 ax=plot(obj.fitData.nonLaserPseudoR2{n},'.:');
            %                 ax.MarkerSize=20;
            %
            %             end
            %             ylim([0 1]); set(gca,'box','off'); ylabel('McFadden pseudo-R^2'); xlabel('Session number');
            %             legend(obj.names(1:end-1));
            
        end
        
        function plotData_psych(obj)
            %Plot summary of all psychomatrices for all mice
            figure('name','mega non-laser plot','color','w');
            numSubjects = length(obj.names)-1; %exclude allSubj
            maxNumSessions = max(cellfun(@length,obj.expRefs(1:numSubjects)));
            for n = 1:numSubjects
                numSessions = max(obj.data{n}.sessionID);
                
                i = 1;
                for s = 1:numSessions
                    D = getrow(obj.data{n},obj.data{n}.sessionID==s & obj.data{n}.laserIdx==0);
                    
                    cVal = unique(unique(D.stimulus(:,1:2)));
                    prop = nan(4,4,3);
                    for cl = 1:length(cVal)
                        for cr = 1:length(cVal)
                            r = D.response(D.stimulus(:,1)==cVal(cl) & D.stimulus(:,2)==cVal(cr));
                            prop(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                        end
                    end
                    
                    for r = 1:3
                        subplot(maxNumSessions,numSubjects*4,(s-1)*(4*numSubjects) + (r) + (n-1)*4);
                        %                         pcolor([cVal;1],[cVal;1],prop(:,:,r)); shading('flat'); caxis([0 1]); axis square;
                        imagesc(prop(:,:,r)); caxis([0 1]); axis square; set(gca,'ydir','normal');
                        set(gca,'box','off','xtick','','ytick','','xcolor','w','ycolor','w');
                        
                        if r==2 && s==1
                            title(obj.names{n});
                        end
                    end
                    
                    %                     pcolor([cVal;1],[cVal;1],prop(:,:,1)); shading('flat'); caxis([0 1]); axis square;
                    %                     set(gca,'box','off','xtick','','ytick','');
                    %                     keyboard;
                end
                
                
                %                 error('incomplete');
                %                 keyboard;
            end
        end
        
        function plotData_Z(obj,varargin)
            %Plot empirical Z space for each
            %contrast condition
            
            if isempty(varargin)
                D = obj.data{end};
                subj = length(obj.names);
            else
                subj = varargin{1};
                D = obj.data{subj};
            end
            
% %             keyboard;
%             [u1,u2,u3]=unique(D.stimulus,'rows');
%             f = tabulate(u3); 
%             f = f(:,3);
%             idx = f>1;
%             %Get only most frequent stimulus types
%             cVal = unique(unique(D.stimulus(idx,:)));
%             
            cVal = unique(unique(D.stimulus));
            
            xlab = 'log(pR/pNG)';
            ylab = 'log(pL/pNG)';
            
            figure('color','w');
            p=nan(length(cVal),length(cVal),3);
            C=nan(length(cVal),length(cVal),2);
            ZL_ci=nan(length(cVal),length(cVal),2);
            ZR_ci=nan(length(cVal),length(cVal),2);
            for cl = 1:length(cVal)
                for cr = 1:length(cVal)
                    C(cl,cr,:) = cVal([cl cr]);
                    r = D.response(D.stimulus(:,1)==cVal(cl) & D.stimulus(:,2)==cVal(cr));
                    p(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                    
%                     psim = squeeze(p(cl,cr,:));
%                     n = length(r);
                    %                         psimDraws=mnrnd(n,psim,5000)/n;
                    %
                    %                         ZL_ci(cl,cr,:) = quantile(log(psimDraws(:,1)./psimDraws(:,3)),[0.025 0.975]);
                    %                         ZR_ci(cl,cr,:) = quantile(log(psimDraws(:,2)./psimDraws(:,3)),[0.025 0.975]);
                    %
                    %                         ZGO_ci(cl,cr,:) = quantile(log((1-psimDraws(:,3))./psimDraws(:,3)),[0.025 0.975]);
                    %                         psimDraws_LR = bsxfun(@rdivide,psimDraws(:,1:2),sum(psimDraws(:,1:2),2));
                    %                         ZLRG_ci(cl,cr,:) = quantile(log(psimDraws_LR(:,1)./psimDraws_LR(:,2)),[0.025 0.975]);
                    %
                    
                end
            end
            
            h1(1)=subplot(2,1,1);
            
            %plot in log L/NG and log R/NG space, irrespective of
            %whether the model was fit in that space or the nested
            %logit space
            
            %                 if obj.nestedModel_flag == 0
            ZL = log(p(:,:,1)./p(:,:,3));
            ZR = log(p(:,:,2)./p(:,:,3));
            %
            %                 ZL_ciLOW = ZL_ci(:,:,1);
            %                 ZL_ciHIGH = ZL_ci(:,:,2);
            %                 ZR_ciLOW = ZR_ci(:,:,1);
            %                 ZR_ciHIGH = ZR_ci(:,:,2);
            %
%             scatter(ZR(:),ZL(:),50,'o'); hold on;
            
%             keyboard;
            plot(ZR(1),ZL(1),'o');
             hold on;
            for cl = 1:length(cVal)
                for cr = 1:length(cVal)
                    cont = cVal([cl cr])/max(cVal);
                    zr = ZR(cl,cr);
                    zl = ZL(cl,cr);
                    
                    if ~isnan(zr) && ~isnan(zl) && ~isinf(zr) && ~isinf(zl)
                        r = fill([zr zr zr+0.1*cont(2) zr ],[zl+0.1 zl-0.1 zl zl+0.1],'r'); r.LineStyle='none';
                        l = fill([zr zr zr-0.1*cont(1) zr ],[zl+0.1 zl-0.1 zl zl+0.1],'b'); l.LineStyle='none';
                    end
                end
            end
            
            for cl = 1:length(cVal)
                idx = ~( isnan(ZL(cl,:)) | isinf(ZL(cl,:)) |  isnan(ZR(cl,:)) | isinf(ZR(cl,:)));
                plot(ZR(cl,idx),ZL(cl,idx),'k:');
            end
            for cr = 1:length(cVal)
                idx = ~( isnan(ZL(:,cr)) | isinf(ZL(:,cr)) |  isnan(ZR(:,cr)) | isinf(ZR(:,cr)));
                plot(ZR(idx,cr),ZL(idx,cr),'k:');
            end
           
            try
                ezplot('2*exp(x)=1+exp(x)+exp(y)');
                ezplot('2*exp(y)=1+exp(x)+exp(y)');
                ezplot('2=1+exp(x)+exp(y)');
            catch
                warning('error on ezplot of contours');
            end
            axis equal;
            xlabel(['empirical ' xlab]); ylabel(['empirical ' ylab]); 
            
            %overlay black dot for random  behaviour
            r = D.response;
            r=sum([r==1 r==2 r==3],1)/length(r);
            hold on;
            h = plot(log(r(2)/r(3)),log(r(1)/r(3)),'k.');
            hold off;
            h.MarkerSize=30;
            
            
            h1(2)=subplot(2,1,2); hold on;
            try
                ezplot('2*exp(x)=1+exp(x)+exp(y)');
                ezplot('2*exp(y)=1+exp(x)+exp(y)');
                ezplot('2=1+exp(x)+exp(y)');
            catch
                warning('error on ezplot of contours');
            end
            param = median(obj.fitData.params{subj},1);
            
            [cr,cl]=meshgrid(cVal);
            cr = cr(:); cl = cl(:);
            
            prop = obj.calculatePhat(obj.model.ZL,obj.model.ZR,param,[cl cr]);
            
            ZL = log(prop(:,1)./prop(:,3)); ZL = reshape(ZL,length(cVal),length(cVal));
            ZR = log(prop(:,2)./prop(:,3)); ZR = reshape(ZR,length(cVal),length(cVal));
            
%             scatter(ZR(:),ZL(:),50,'o');
            plot(ZR(1),ZL(1),'o');
            plot(ZR,ZL,'k:');
            plot(ZR',ZL','k:');

            for cl = 1:length(cVal)
                for cr = 1:length(cVal)
                    cont = cVal([cl cr])/max(cVal);
                    zr = ZR(cl,cr);
                    zl = ZL(cl,cr);
                    
                    if ~isnan(zr) && ~isnan(zl) && ~isinf(zr) && ~isinf(zl)
                        r = fill([zr zr zr+0.1*cont(2) zr ],[zl+0.1 zl-0.1 zl zl+0.1],'r'); r.LineStyle='none';
                        l = fill([zr zr zr-0.1*cont(1) zr ],[zl+0.1 zl-0.1 zl zl+0.1],'b'); l.LineStyle='none';
                    end
                end
            end
            
            
            axis equal;
            xlabel(['predicted ' xlab]); ylabel(['predicted ' ylab]); title('');
            
            linkaxes(h1,'xy')
        end
    end
    
    methods (Access=private)
        
        function eRefs = getExpRefs(~,name,expt)
            try
                load(['\\basket.cortexlab.net\home\glmmulti_files\' expt '_' name '_EXPREFS.mat']);
            catch
                warning('Data not found, recalculating now...');
                
                eRefs = vertcat(dat.listExps(name));
                good = ones(length(eRefs),1);
                
                for b = 1:length(eRefs)
                    
                    try
                        p = load(dat.expFilePath(eRefs{b},'parameters','m'));
                        g = GLM(eRefs{b}); %check if any lasermanip problems
                        if exist(dat.expFilePath(eRefs{b},'lasermanip','m'), 'file') == 2
                            l = laserGLM(eRef);
                        end
                        
                        %remove any sessions with performance < 75%
                        if mean(g.data.feedbackType==1)<0.75
                            good(b)=0;
                        end
                        
                        switch(expt)
                            case '2D'
                                stim = p.parameters.visCueContrast;
                                if ~any(min(stim,[],1)>0)
                                    good(b)=0;
                                end
                                
                                if length(g.data.response)<150
                                    good(b)=0;
                                end
                                
                            case '1D'
                                stim = p.parameters.visCueContrast;
                                if any(min(stim,[],1)>0)
                                    good(b)=0;
                                end
                                
                                if length(l.data.response)<150
                                    good(b)=0;
                                end
                                
                            otherwise
                                error('choose something');
                        end
                        
                        if good(b)==1
                            eRefs{b}
                        end
                        
                    catch
                        good(b)=0;
                    end
                end
                
                eRefs = eRefs(good==1);
                
                save(['\\basket.cortexlab.net\home\glmmulti_files\' expt '_' name '_EXPREFS.mat'],'eRefs');
            end
        end
        
        function [p,exitflag] = util_optim(obj,ZL,ZR,cont,resp,LB,UB)
            objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL,ZR,cont,resp) + obj.fitLambda*sum(abs(PARAMETERS).^2) );
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000,'Display','off');
            [p,~,exitflag] = fmincon(objective,zeros(1,length(LB)), [], [], [], [], LB, UB, [], options);

            if exitflag<1
                warning('Did not converge, trying Opti');
                options = optiset('display','final','solver','NOMAD');
                Opt = opti('fun',objective,'x0',zeros(1,length(LB)),'bounds',LB,UB,'options',options);
                [pTemp,~,exitflag] = Opt.solve;
                p = pTemp';
            end
            
            if exitflag<1
                warning('Failed to converge with opti. Fix this!');
                keyboard;
            end
        end
        
        function logLik = calculateLogLik(obj,testParams,ZL,ZR,inputs,responses)
            phat = obj.calculatePhat(ZL,ZR,testParams,inputs);
            p = (responses==1).*phat(:,1) + (responses==2).*phat(:,2) + (responses==3).*phat(:,3);
            p(p<0)=0; 
            logLik = nanmean(log2(p));

        end
        
        function phat = calculatePhat(obj,ZL,ZR,testParams,inputs)
            zl = ZL(testParams,inputs);
            zr = ZR(testParams,inputs);
            pL = exp(zl)./(1+exp(zl)+exp(zr));
            pR = exp(zr)./(1+exp(zl)+exp(zr));
            pNG = 1 - pL - pR;
            
            
            if obj.nestedModel_flag==1 %then ZL and ZR were actually ZGO/NG and ZL/R
                zgo = zl;
                zlr = zr;
                
                pGO = exp(zgo)./(1+exp(zgo));
                pL = exp(zlr)./(1+exp(zlr));
                pR = 1-pL;
                pL = pL.*pGO; pR = pR.*pGO;
                pNG = 1-pGO;
            end
            
            phat = [pL pR pNG];
            
        end

    end
end