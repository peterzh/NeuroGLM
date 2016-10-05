classdef omnibusLaserGLM
    properties
        expt;
        names;
        expRefs;
        inactivationCoords;
        data;
        fitData;
        fitMethod='fmincon';
        guess_bpt = [];
    end
    
    properties (Access=private)
        cfn_parameters;
        significance_flag = 0;
        bilateral_flag = 0;
        moreData = 0;
    end
    
    methods
        function obj = omnibusLaserGLM(expt,nameList)
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
            obj.inactivationCoords = cell(1,numSubjects);
            obj.data = cell(1,numSubjects);
            obj.cfn_parameters = cell(1,numSubjects);
            for n = 1:numSubjects
                D = struct;
                c50_parameters = [];

                for session = 1:length(obj.expRefs{n})
                    eRef = obj.expRefs{n}{session};
                    try
                        L=load(dat.expFilePath(eRef, 'laserManip', 'm'));
                        L.laserCoordByTrial;
                    catch
                        getLaserLabels(inputData);
                        L=load(dat.expFilePath(inputData, 'laserManip', 'm'));
                        L.laserCoordByTrial; %Test if this variable exists, that's all...
                    end
                    
                    %if bilateral then replace each row by its virtual
                    %coordinate
                    try
                        L.coordList_unadjusted(1,:); %in bilat expts this field exists
                        for i=1:size(L.laserCoordByTrial,1)
                            if ~isnan(L.laserCoordByTrial(i,1))
                                matchIdx = sum(abs(bsxfun(@minus,L.laserCoordByTrial(i,:),L.coordList)),2)==0;
                                L.laserCoordByTrial(i,:) = [L.coordList_unadjusted(matchIdx,:) 0];
                            end
                        end
                        
                    catch
                    end
                    
                    if size(L.laserCoordByTrial,2)==2
                        L.laserCoordByTrial = [L.laserCoordByTrial zeros(size(L.laserCoordByTrial,1),1)];
                    end
                    
                    d = struct('contrast_cond',[],'response',[],'repeatNum',[],'feedbackType',[],'RT',[]);
                    block = dat.loadBlock(eRef);
                    trials = block.trial;
                    
                    for t=1:block.numCompletedTrials
                        d.contrast_cond(t,:) = trials(t).condition.visCueContrast';
                        d.response(t,1) = trials(t).responseMadeID';
                        d.repeatNum(t,1) = trials(t).condition.repeatNum;
                        d.feedbackType(t,1) = trials(t).feedbackType;
                        d.RT(t,1) = trials(t).responseMadeTime-trials(t).interactiveStartedTime;
                    end
                    
                    d.laserCoord = L.laserCoordByTrial;

                    if obj.moreData == 1
                        pupilenergy = obj.addPupilData(n,obj.expt);
                        d.pupilenergy = pupilenergy{session};
                    end
                    d = structfun(@(x)(x(6:(end-14),:)),d,'uni',0); %trim first 5 trials and last 15
                    
%                     %Fit simple C50 models to estimate c50 & n params
                    e = getrow(d,isnan(d.laserCoord(:,1)));
                    g=GLM(e).setModel('C50-subset').fit;
                    c50_parameters(session,:) = g.parameterFits(5:6);
%                     
                    d.sessionID = ones(length(d.response),1)*session;
                    D = addstruct(D,d);
                end
                
                D = getrow(D,D.repeatNum==1);
                
                %create laserIdx variable
                [sites,~,ic] = unique(D.laserCoord,'rows');
                nanIdx = find(isnan(sites(:,1)),1,'first');
                inactivationSite = sites(1:nanIdx-1,:);
                ic(ic>=nanIdx)=0;
                D.laserIdx = ic;

                obj.inactivationCoords{n} = inactivationSite;
                obj.data{n} = D;
                
                tab = tabulate(D.response);
                tab = tab(:,3)/100;
                obj.guess_bpt(n)=sum(tab.*log2(tab));
                
                obj.cfn_parameters{n} = mean(c50_parameters,1);
            end
            
            if ~isempty(strfind(expt,'bilateral'))
                obj.bilateral_flag = 1;
            end
   
            
        end
        
        function [ZL,ZR,numP,offsetMode] = getModel(obj,model,biasMode,subj)
            cfn = @(c)(obj.cfn(c,subj));
            switch(biasMode)
                case 'biasAsOffset'
                    switch(model)
                        case 'offsetOnly'
                            ZL = @(off,p,c)(off(:,1) );
                            ZR = @(off,p,c)(off(:,2) );
                            numP = NaN;
                        case 'bias'
                            ZL = @(off,p,c)(off(:,1) + p(1) );
                            ZR = @(off,p,c)(off(:,2) + p(2) );
                            numP = 2;
                        case 'sub_sens'
                            ZL = @(off,p,c)(off(:,1) + p(1).*cfn(c(:,1)) );
                            ZR = @(off,p,c)(off(:,2) + p(2).*cfn(c(:,2)) );
                            numP = 2;
                        case 'bias+sub_sens'
                            ZL = @(off,p,c)(off(:,1) + p(1) + p(3).*cfn(c(:,1)) );
                            ZR = @(off,p,c)(off(:,2) + p(2) + p(4).*cfn(c(:,2)) );
                            numP = 4;
                    end
                    offsetMode = 'Z';
                    
                case 'biasAsContrast'
                    switch(model)
                        case 'offsetOnly'
                            ZL = @(p_offset,p,c)( p_offset(:,3).*(cfn(c(:,1)) + p_offset(:,1)) );
                            ZR = @(p_offset,p,c)( p_offset(:,4).*(cfn(c(:,2)) + p_offset(:,2)) );
                            numP = NaN;
                        case 'bias+sub_sens'
                            ZL = @(p_offset,p,c)( (p_offset(:,3) + p(3)).*(cfn(c(:,1)) + p_offset(:,1) + p(1)) );
                            ZR = @(p_offset,p,c)( (p_offset(:,4) + p(4)).*(cfn(c(:,2)) + p_offset(:,2) + p(2)) );
                            numP = 4;
                        case 'bias'
                            ZL = @(p_offset,p,c)( ( p_offset(:,3) ).*(cfn(c(:,1)) + p_offset(:,1) + p(1)) );
                            ZR = @(p_offset,p,c)( ( p_offset(:,4) ).*(cfn(c(:,2)) + p_offset(:,2) + p(2)) );
                            numP = 2;
                        case 'sub_sens'
                            ZL = @(p_offset,p,c)( ( p_offset(:,3) + p(1) ).*(cfn(c(:,1)) + p_offset(:,1) ) );
                            ZR = @(p_offset,p,c)( ( p_offset(:,4) + p(2) ).*(cfn(c(:,2)) + p_offset(:,2) ) );
                            numP = 2;
                    end
                    offsetMode = 'P';
                    
            end
        end
        
        function obj = fit(obj,lambda)
            f=figure;
            obj.fitData.nonLaserModel = 'bias+sub_sens';
            obj.fitData.laserModel = 'bias+sub_sens';
            obj.fitData.biasMode = 'biasAsContrast';
            obj.fitData.lambda = lambda;
            
            numSubjects = length(obj.names);
            %             numSubjects = 1;
            obj.fitData.params = cell(1,numSubjects);
            obj.fitData.nonLaserParams = cell(1,numSubjects);
            for n = 1:numSubjects
                
                %First fit non-laser portion of the data
                
                %Manual optim method per-session
                [ZL_nL,ZR_nL,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                nL = obj.data{n}.laserIdx==0;
                
                %Define which parameters are to be regularised
                if numP == 4
                                        reg = [3:4];
%                     reg = 1:numP;
                else
                    reg = 1:numP;
                end
                
                sessions = unique(obj.data{n}.sessionID(nL));
                p_nL = nan(length(sessions),numP);
                for s = 1:length(sessions)
                    cont = obj.data{n}.contrast_cond(nL & obj.data{n}.sessionID==sessions(s),:);
                    resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                    
                    offset = zeros(length(resp),4);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS(reg).^2) );
                    
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                            p_nL(s,:) = Opt.solve';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                            p_nL(s,:) = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                    end
                end
                
                obj.fitData.nonLaserParams{n} = p_nL;
                
                %Then go through each inactivation location and fit
                %extended models
                p_L = []; loglik=[];
                for site = 1:size(obj.inactivationCoords{n},1)
                    L_idx = obj.data{n}.laserIdx==site;
                    cont = obj.data{n}.contrast_cond(L_idx,:);
                    resp = obj.data{n}.response(L_idx);
                    sess = obj.data{n}.sessionID(L_idx);
                    
                    %Calculate offset from non-laser fit
                    %Offset from optim method
                    
                    switch (offsetMode)
                        case 'P'
                            %offset on parameters
                            offset = p_nL(sess,:);
                        case 'Z'
                            %offset on ZL and ZR
                            offset = [];
                            for t = 1:sum(L_idx)
                                offset(t,:) = [ZL_nL([0 0],p_nL(sess(t),:),cont(t,:)),...
                                    ZR_nL([0 0],p_nL(sess(t),:),cont(t,:)) ];
                            end
                    end
                    
                    %Use that offset to manually fit a model at each
                    %inactivation site
                    [ZL,ZR,numP,offsetMode] = obj.getModel(obj.fitData.laserModel,obj.fitData.biasMode,n);
                    objective = @(PARAMETERS) ( -obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS(reg).^2) );
                    
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                            [p,loglik(site,1),~,~] = Opt.solve;
                            p_L(site,:) = p';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                            [p_L(site,:),loglik(site,1),~,~] = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                    end
                    
                    %                     obj.fitData.siteX{site,n} = X;
                    %                     obj.fitData.siteFit{site,n} = fit;
                end
                obj.fitData.params{n} = p_L;
                %                 obj.fitData.paramsLL{n} = -loglik;
            end
            
            disp('Done!'); close(f);
        end
        
        function obj = nonLaserBootstrap(obj)
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
                    cont = obj.data{n}.contrast_cond(nL_bootstrap,:);
                    resp = obj.data{n}.response(nL_bootstrap);
                    
                    offset = zeros(length(resp),4);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS.^2) );
                    
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
        %
        function obj = fitNull(obj)
            figure;
            kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7); hold on;
            numSites = size(obj.inactivationCoords{1},1);
            cols = [linspace(0,1,numSites)' linspace(1,0,numSites)' zeros(numSites,1) ];
            cols(obj.inactivationCoords{1}(:,2)>0,3) = 1;
            scatter(obj.inactivationCoords{1}(:,2),obj.inactivationCoords{1}(:,1),200,cols,'filled');
            
            numSubjects = length(obj.names);
            numIter = 500;
            
            obj.fitData.paramsNull = cell(1,numSubjects);
            figure('name','Sampling from null distribution of laser parameters');
            for n = 1:numSubjects
                %Get all non-laser trials
                tab = tabulate(obj.data{n}.laserIdx);
                proportion = mean(tab(2:end,3)/100);
                
                subplot(1,numSubjects,n); title(obj.names{n}); hold on;
                
                numLP = size(obj.fitData.params{n},2);
                if numLP == 4
                    scatter(obj.fitData.params{n}(:,1),obj.fitData.params{n}(:,3),30,cols,'filled');
                elseif numLP == 2
                    scatter(obj.fitData.params{n}(:,1),obj.fitData.params{n}(:,2),30,cols,'filled');
                end
                
                loglik = [];
                for iter = 1:numIter
                    disp(iter);
                    nL = obj.data{n}.laserIdx==0;
                    %In each interation, assign a proportion of them to be
                    %a laser site. Ensure it is a stratified sample of choice!
                    numInSample = round(proportion*sum(nL));
                    strat_response = tabulate(obj.data{n}.response(nL));
                    
                    q_idx = [];
                    for r = 1:3
                        numOfChoice=round(numInSample*strat_response(r,3)/100);
                        IDxOfChoice=(obj.data{n}.response==r & nL);
                        q_idx = [q_idx; randsample(find(IDxOfChoice),numOfChoice)];
                    end
                    nL(q_idx) = 0; %Set these non-laser trials to be a laser trial
                    
                    %Of the remaining nonLaser trials, fit the nonLaser mod
                    
                    %Manual optim method per-session
                    [ZL_nL,ZR_nL,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                    %Define which parameters are to be regularised
                    if numP == 4
                        reg = [3:4];
                    else
                        reg = 1:numP;
                    end
                    
                    sessions = unique(obj.data{n}.sessionID(nL));
                    p_nL = nan(length(sessions),numP);
                    for s = 1:length(sessions)
                        cont = obj.data{n}.contrast_cond(nL & obj.data{n}.sessionID==sessions(s),:);
                        resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                        
                        offset = zeros(length(resp),4);
                        objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS(reg).^2));
                        
                        switch(obj.fitMethod)
                            case 'opti'
                                options = optiset('display','final','solver','NOMAD');
                                Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                                p_nL(s,:) = Opt.solve';
                                
                            case 'fmincon'
                                options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000,'Display','off');
                                p_nL(s,:) = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                        end
                    end
                    
                    cont = obj.data{n}.contrast_cond(q_idx,:);
                    resp = obj.data{n}.response(q_idx);
                    sess = obj.data{n}.sessionID(q_idx);
                    
                    switch (offsetMode)
                        case 'P'
                            %offset on parameters
                            offset = p_nL(sess,:);
                        case 'Z'
                            %offset on ZL and ZR
                            offset = [];
                            for t = 1:length(q_idx)
                                offset(t,:) = [ZL_nL(0,p_nL(sess(t),:),cont(t,:)),...
                                    ZR_nL(0,p_nL(sess(t),:),cont(t,:)), 0, 0 ];
                            end
                    end
                    
                    [ZL,ZR,numP,offsetMode] = obj.getModel(obj.fitData.laserModel,obj.fitData.biasMode,n);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS(reg).^2));
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                            [p,loglik(iter,1)] = Opt.solve;
                            obj.fitData.paramsNull{n}(iter,:) = p';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000,'Display','off');
                            [p,loglik(iter,1)] = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                            obj.fitData.paramsNull{n}(iter,:) = p;
                    end
                    if numLP == 4
                        plot(p(1),p(3),'k.','markersize',10); drawnow; xlabel('B_L'); ylabel('S_L');
                    elseif numLP == 2
                        plot(p(1),p(2),'k.','markersize',10); drawnow; xlabel('B_L'); ylabel('B_R');
                    end
                    
                end
                
                %                 keyboard;
                
                obj.fitData.paramsNullLL{n}=-loglik;
                
                %                 figure; gplotmatrix(obj.fitData.paramsNull{n});
            end
            
            obj.significance_flag = 1;
        end
        
        function obj = clearNull(obj)
            obj.significance_flag=0;
            try
                obj.fitData = rmfield(obj.fitData,{'paramsNull','paramsNullLL'});
            catch
                warning('Null data does not exist');
            end
        end
        
        function plotFit_Null(obj)
            figure;
            kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7); hold on;
            numSites = size(obj.inactivationCoords{1},1);
            cols = [linspace(0,1,numSites)' linspace(1,0,numSites)' zeros(numSites,1) ];
            cols(obj.inactivationCoords{1}(:,2)>0,3) = 1;
            %             cols(:,3) = obj.inactivationCoords{1}(:,2)/7 + .5;
            scatter(obj.inactivationCoords{1}(:,2),obj.inactivationCoords{1}(:,1),200,cols,'filled');
            %             text(obj.inactivationCoords{1}(:,2),obj.inactivationCoords{1}(:,1),cellfun(@(c)num2str(c),num2cell(1:numSites),'uni',0));
            
            numLP = size(obj.fitData.paramsNull{1},2);
            
            figure('name','Sampling from null distribution of laser parameters','color','w');
            numSubjects = length(obj.names);
            for n = 1:numSubjects
                
                if numLP == 4
                    %Plot original B vs S
                    subplot(2,numSubjects,n); title(obj.names{n}); hold on;
                    plot(obj.fitData.paramsNull{n}(:,1),obj.fitData.paramsNull{n}(:,3),'.','markerfacecolor',[1 1 1]*0.5,'markeredgecolor',[1 1 1]*0.5);
                    scatter(obj.fitData.params{n}(:,1),obj.fitData.params{n}(:,3),30,cols,'filled');
                    xlabel('\Delta B_L'); ylabel('\Delta S_L'); axis square;
                    %                     t=text(obj.fitData.params{n}(:,1) + 0.05,obj.fitData.params{n}(:,3),cellfun(@(c)num2str(c),num2cell(1:numSites),'uni',0));
                    %                     set(t,'FontSize',7);
                    
                    subplot(2,numSubjects,numSubjects+n); hold on;
                    plot(obj.fitData.paramsNull{n}(:,2),obj.fitData.paramsNull{n}(:,4),'.','markerfacecolor',[1 1 1]*0.5,'markeredgecolor',[1 1 1]*0.5);
                    scatter(obj.fitData.params{n}(:,2),obj.fitData.params{n}(:,4),30,cols,'filled');
                    xlabel('\Delta B_R'); ylabel('\Delta S_R'); axis square;
                    %                     t=text(obj.fitData.params{n}(:,2) + 0.05,obj.fitData.params{n}(:,4),cellfun(@(c)num2str(c),num2cell(1:numSites),'uni',0));
                    %                     set(t,'FontSize',7);
                    
                elseif numLP == 2
                    subplot(1,numSubjects,n); title(obj.names{n}); hold on;
                    plot(obj.fitData.paramsNull{n}(:,1),obj.fitData.paramsNull{n}(:,2),'.','markerfacecolor',[1 1 1]*0.5,'markeredgecolor',[1 1 1]*0.5);
                    scatter(obj.fitData.params{n}(:,1),obj.fitData.params{n}(:,2),30,cols,'filled');
                    xlabel('B_L'); ylabel('B_R'); axis square;
                    
                    %                     t=text(obj.fitData.params{n}(:,1) + 0.05,obj.fitData.params{n}(:,2),cellfun(@(c)num2str(c),num2cell(1:numSites),'uni',0));
                    %                     set(t,'FontSize',7);
                end
            end
            
            set(get(gcf,'children'),'xlim',[-1 1]*0.5,'ylim',[-1 1]*10);
        end
        
        function crossval_LaserEffect(obj)
            numFolds = 10;
            laserModels = {'bias','bias+sub_sens'};
            numSubjects = length(obj.names);
            
            figure('name',[num2str(numFolds) '-fold CV ' obj.fitData.biasMode],'color','w');
            
            logLik = cell(1,numSubjects);
            for n = 1:numSubjects
                cv = cvpartition(obj.data{n}.response,'KFold',numFolds);
                phat = nan(length(obj.data{n}.response),length(laserModels)+1);
                for fold = 1:cv.NumTestSets
                    
                    trainC = obj.data{n}.contrast_cond(cv.training(fold),:);
                    trainR = obj.data{n}.response(cv.training(fold));
                    trainL = obj.data{n}.laserIdx(cv.training(fold));
                    trainSID = obj.data{n}.sessionID(cv.training(fold));
                    
                    %Fit non-laser model to training set
                    [ZL_nL,ZR_nL,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                    nL = (trainL==0);
                    
                    sessions = unique(trainSID(nL));
                    p_nL = nan(length(sessions),numP);
                    disp('NON-LASER TRAIN LOGLIK RELATIVE TO GUESS');
                    for s = 1:length(sessions)
                        cont = trainC(nL & trainSID==sessions(s),:);
                        resp = trainR(nL & trainSID==sessions(s));
                        
                        offset = zeros(length(resp),4);
                        objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS.^2) );
                        
                        switch(obj.fitMethod)
                            case 'opti'
                                options = optiset('display','final','solver','NOMAD');
                                Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                                p_nL(s,:) = Opt.solve';
                                
                            case 'fmincon'
                                options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000,'Display','off');
                                p_nL(s,:) = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                        end
                        
                        %debug: check fits were of good quality
                        p=obj.calculatePhat(ZL_nL,ZR_nL,offset,p_nL(s,:),cont);
                        loglik=mean(log2(p(:,1).*(resp==1) + p(:,2).*(resp==2) + p(:,3).*(resp==3)))-obj.guess_bpt(n);
                        disp(['Session ' num2str(s) ': ' num2str(loglik)]);
                    end
                    
                    %Now fit laser part to the training set
                    p_L = cell(size(obj.inactivationCoords{n},1),length(laserModels));
                    disp('LASER TRAIN LOGLIK RELATIVE TO GUESS');
                    for site = 1:size(obj.inactivationCoords{n},1)
                        disp(['Site ' num2str(site)]);
                        cont = trainC(trainL==site,:);
                        resp = trainR(trainL==site);
                        sid = trainSID(trainL==site);
                        
                        switch(offsetMode)
                            case 'P'
                                offset = p_nL(sid,:);
                            case 'Z'
                                [ZL,ZR,~,~] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                                offset = [];
                                for t = 1:sum(trainL==site)
                                    offset(t,:) = [ZL([0 0],p_nL(sid(t),:),cont(t,:)),...
                                        ZR([0 0],p_nL(sid(t),:),cont(t,:)),0,0];
                                end
                        end
                        [ZL,ZR,~,~] = obj.getModel('offsetOnly',obj.fitData.biasMode,n);
                        
                        p = obj.calculatePhat(ZL,ZR,offset,[],cont);
                        loglik=mean(log2(p(:,1).*(resp==1) + p(:,2).*(resp==2) + p(:,3).*(resp==3)))-obj.guess_bpt(n);
                        disp(['    non-laser model : ' num2str(loglik)]);
                        
                        for model = 1:length(laserModels)
                            [ZL,ZR,numP,offsetMode] = obj.getModel(laserModels{model},obj.fitData.biasMode,n);
                            objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp));
                            
                            switch(obj.fitMethod)
                                case 'opti'
                                    options = optiset('display','final','solver','NOMAD');
                                    Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                                    p_L{site,model} = Opt.solve';
                                    
                                case 'fmincon'
                                    options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000,'Display','off');
                                    p_L{site,model} = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                            end
                            
                            %debug: check fits were of good quality
                            p=obj.calculatePhat(ZL,ZR,offset,p_L{site,model},cont);
                            loglik=mean(log2(p(:,1).*(resp==1) + p(:,2).*(resp==2) + p(:,3).*(resp==3)))-obj.guess_bpt(n);
                            disp(['    ' laserModels{model} ' : ' num2str(loglik)]);
                        end
                    end
                    
                    %Now test these fits at the sites in the test set
                    testIdx = find(cv.test(fold));
                    testC = obj.data{n}.contrast_cond(testIdx,:);
                    testR = obj.data{n}.response(testIdx);
                    testL = obj.data{n}.laserIdx(testIdx);
                    testSID = obj.data{n}.sessionID(testIdx);
                    
                    %go through each trial in the test set
                    for t = 1:length(testR)
                        cont = testC(t,:);
                        resp = testR(t);
                        las = testL(t);
                        sid = testSID(t);
                        
                        switch(offsetMode)
                            case 'P'
                                offset = p_nL(sid,:);
                            case 'Z'
                                [ZL,ZR,~,~] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                                offset = [ZL([0 0],p_nL(sid,:),cont),...
                                    ZR([0 0],p_nL(sid,:),cont),0,0];
                        end
                        [ZL,ZR,~,~] = obj.getModel('offsetOnly',obj.fitData.biasMode,n);
                        
                        
                        %Test the prediction of each laser model on that
                        %trial (if a laser trial)
                        if las>0
                            %non-laser model
                            
                            p = obj.calculatePhat(ZL,ZR,offset,[],cont);
                            phat(testIdx(t),1) = (resp==1)*p(1) + (resp==2)*p(2) + (resp==3)*p(3);
                            
                            %other laser models
                            for model = 1:length(laserModels)
                                [ZL,ZR,~,offsetMode] = obj.getModel(laserModels{model},obj.fitData.biasMode,n);
                                p = obj.calculatePhat(ZL,ZR,offset,p_L{las,model},cont);
                                phat(testIdx(t),model+1) = (resp==1)*p(1) + (resp==2)*p(2) + (resp==3)*p(3);
                            end
                        end
                    end
                end
                
                subplot(2,length(obj.names),n);
                bar(nanmean(log2(phat))-obj.guess_bpt(n));
                set(gca,'XTickLabel',['noLaser' laserModels],'XTickLabelRotation',90,'box','off');
                title(obj.names{n}); ylabel('CV loglik relative to guess');
                
                subplot(2,length(obj.names),n+length(obj.names));
                for m = 1:size(phat,2)
                    logLik{n}(:,m) = pivottable(obj.data{n}.laserIdx,[],log2(phat(:,m)),'nanmean');
                end
                logLik{n}(1,:) = [];
                
                %identify which was the winner among the models
                [~,winner] = max(logLik{n},[],2);
                
                delta_LL = logLik{n}(:,3) - logLik{n}(:,1);
                %                 delta_LL(winner==1) = NaN;
                
                ap = obj.inactivationCoords{n}(:,1);
                ml = obj.inactivationCoords{n}(:,2);
                if obj.bilateral_flag==1
                    ap = [ap; ap];
                    ml = [ml; -ml];
                    %                     winner = [winner;winner];
                    delta_LL = [delta_LL;delta_LL];
                end
                
                scatter(ml,ap,200,delta_LL,'s','filled'); caxis([-1 1]*max(abs(delta_LL)));
                axis equal; set(gca,'xtick','','ytick','','xcolor','w','ycolor','w');
                cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                colormap(flipud(cmap)); colorbar; title('LL(Bias+Sens) - LL(NoLaser)');
                drawnow;
            end
            %
            % %             keyboard;
            %             %Plot average model winner across mice
            %             logLik_ave = mean(cat(4,logLik{:}),4);
            %             [~,winner] = max(logLik_ave,[],2);
            %             figure('name','Average','color','w');
            %
            %             ap = obj.inactivationCoords{1}(:,1);
            %             ml = obj.inactivationCoords{1}(:,2);
            %             if obj.bilateral_flag==1
            %                 ap = [ap; ap];
            %                 ml = [ml; -ml];
            %                 winner = [winner;winner];
            %             end
            %
            %             scatter(ml,ap,200,cols(winner,:),'s','filled'); axis equal; set(gca,'xtick','','ytick','','xcolor','w','ycolor','w');
        end
        
        function plotFit_NonLaser(obj)
            numSubjects = length(obj.names);
            
            %Plot session-by-session bias and sensitivity parameters
            figure('name','Session by session non-laser bias and sensitivity changes','color','w');
            for n = 1:numSubjects
                subplot(numSubjects,1,n);
                p_nL = obj.fitData.nonLaserParams{n};
                
                %                 if exist('obj.fitData.nonLaserParamsSTD')
                try
                    errorbar(p_nL,obj.fitData.nonLaserParamsSTD{n});
                catch
                    plot(p_nL,'o-');
                end
                %                 else
                %                     plot(p_nL);
                %                 end
                set(gca,'box','off');
                title(obj.names{n});
                if n == numSubjects
                    xlabel('Session number');
                    legend('B_L','B_R','S_L','S_R');
                end
                %                 keyboard;
            end
            
            %Plot non-laser predictions against actual psychometric data.
            cols = {[0 0.4470 0.7410],...
                [0.8500 0.3250 0.0980],...
                [0.9290 0.6940 0.1250]};
            
            for n = 1:numSubjects
                figure('name',[obj.names{n} ' non-laser predictions against actual psychometric data']);
                [ZL_nL,ZR_nL,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                
                numSessions = max(obj.data{n}.sessionID);
                for s = 1:numSessions
                    p_nL = obj.fitData.nonLaserParams{n}(s,:);
                    idx = obj.data{n}.laserIdx==0 & obj.data{n}.sessionID==s;
                    
                    cont = obj.data{n}.contrast_cond(idx,:);
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
                        ped_c_diff = diff(cont(ped_idx,:),[],2);
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
                        p_hat = obj.calculatePhat(ZL_nL,ZR_nL,zeros(length(testCont),4),p_nL,testCont);
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
                end
                set(gcf,'color','w');
            end
        end
        
        function plotFit_Laser(obj)
            numSubjects = length(obj.names);
            
            %Plot laser parameter maps
            figure('name','Laser parameter maps','color','w');
            for n = 1:numSubjects
                p_L = obj.fitData.params{n};
                
                if size(p_L,2) == 2
                    plotVal = {p_L(:,1),p_L(:,2)};
                    param_labels = {'b_L','b_R'};
                elseif size(p_L,2) == 4
                    plotVal = {p_L(:,1),p_L(:,2),p_L(:,3),p_L(:,4)};
                    param_labels = {'b_L','b_R','s_L','s_R'};
                end
                
                if obj.significance_flag == 1
                    pvalue = obj.sig{n};
                    dotSize = 35 + (pvalue<0.05)*80 + (pvalue<0.01)*100;
                else
                    dotSize = ones(size(p_L))*200;
                end
                
                for i = 1:length(plotVal)
                    subplot(numSubjects,length(plotVal),length(plotVal)*n-length(plotVal) + i);
                    
                    ap = obj.inactivationCoords{n}(:,1);
                    ml = obj.inactivationCoords{n}(:,2);
                    dot = dotSize(:,i);
                    if obj.bilateral_flag==1
                        ap = [ap; ap];
                        ml = [ml; -ml];
                        plotVal{i} = [plotVal{i}; plotVal{i}];
                        dot = [dot;dot];
                    end
                    
                    scatter(ml,ap,dot,plotVal{i},'s','filled'); axis equal;
                    %                     if obj.significance_flag == 1
                    %                         pvalue = obj.sig{n};
                    %                         dotSize = 35 + (pvalue<0.05)*80 + (pvalue<0.01)*100;
                    %                         scatter(ml,ap,dotSize,plotVal{i},'s','filled'); axis equal;
                    %                     else
                    %                         scatter(ml,ap,200,plotVal{i},'s','filled'); axis equal;
                    %                     end
                    
                    caxis([-1 1]*max(abs(plotVal{i}))); colorbar;
                    set(gca,'box','off','ytick','','xtick','','xcolor','w','ycolor','w');
                    
                    if i==1
                        set(gca,'ycolor','k');
                        ylabel(obj.names{n});
                    end
                    if n==1
                        title(param_labels{i});
                    end
                end
                
                cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                colormap(flipud(cmap));
            end
            
            %Plot average laser parameter maps
            figure('name','AVERAGE Laser parameter maps','color','w');
            %             N = cellfun(@(c)(length(c.response)),obj.data); %num trials
            %             weight(1,1,:) = N/sum(N);
            %             weightedparams = bsxfun(@times,weight,cat(3,obj.fitData.params{:}));
            %             p_L = sum(weightedparams,3); %weighted average
            p_L = mean(cat(3,obj.fitData.params{:}),3);
            if size(p_L,2) == 2
                plotVal = {p_L(:,1),p_L(:,2),p_L(:,1)-p_L(:,2), p_L(:,1)+p_L(:,2)};
                param_labels = {'b_L','b_R','b_L-b_R','b_L+b_R'};
            elseif size(p_L,2) == 4
                plotVal = {p_L(:,1),p_L(:,2),p_L(:,3),p_L(:,4)};
                param_labels = {'b_L','b_R','s_L','s_R'};
            end
            
            for i = 1:length(plotVal)
                subplot(1,length(plotVal),i);
                
                ap = obj.inactivationCoords{1}(:,1);
                ml = obj.inactivationCoords{1}(:,2);
                if obj.bilateral_flag==1
                    ap = [ap; ap];
                    ml = [ml; -ml];
                    plotVal{i} = [plotVal{i}; plotVal{i}];
                end
                
                scatter(ml,ap,200,plotVal{i},'s','filled'); axis equal;
                caxis([-1 1]*max(abs(plotVal{i}))); colorbar;
                set(gca,'box','off','ytick','','xtick','','xcolor','w','ycolor','w');
                title(param_labels{i});
            end
            
            cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
            colormap(flipud(cmap));
            
            %Plot psych curves showing the laser effect at a few sites, for
            %each mouse separately
            
            if obj.bilateral_flag==0
                siteLabels = {'Left V1','Right V1','Left M2','Right M2','Left Barrel','Right Barrel'};
                siteID = [10 15 48 49 25 32];
            elseif obj.bilateral_flag==1
                siteLabels = {'V1','M2','Barrel'};
                siteID = [7 24 16];
            end
            
            for n = 1:numSubjects
                figure('name',[obj.names{n} ' laser psych effects']);
                [ZL_nL,ZR_nL,~,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                
                for site = 1:length(siteID)
                    %Plot non-laser curve, use SESSION AVERAGE of parameters for now
                    p_nL = mean(obj.fitData.nonLaserParams{n});
                    p_L = obj.fitData.params{n}(siteID(site),:);
                    
                    cont = obj.data{n}.contrast_cond;
                    cVals = unique(cont(:));
                    numPedestals = length(cVals)-1;
                    
                    for ped = 1:numPedestals
                        subplot(numPedestals,length(siteID),(ped-1)*length(siteID) + site); hold on;
                        
                        %Plot non-laser predictions
                        uC = unique(diff(cont(min(cont,[],2)==cVals(ped),:),[],2));
                        testCont = [linspace(max(abs(uC))+0.1,0,100)' zeros(100,1); zeros(100,1) linspace(0,max(abs(uC))+0.1,100)'] + cVals(ped);
                        p_hat = obj.calculatePhat(ZL_nL,ZR_nL,zeros(length(testCont),4),p_nL,testCont);
                        set(gca,'ColorOrderIndex',1);
                        plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
                        
                        %Try plotting NULL LASER psych bounds
                        try
                            p_null = obj.fitData.paramsNull{n};
                            if size(p_null,2)==2
                                p_null = [p_null, zeros(size(p_null))];
                            end
                            
                            p_totalNull = bsxfun(@plus,p_nL,p_null);
                            numIter = size(p_totalNull,1);
                            p_hat = nan(size(testCont,1),3,numIter);
                            for i = 1:numIter
                                p_hat(:,:,i) = obj.calculatePhat(ZL_nL,ZR_nL,zeros(length(testCont),4),p_totalNull(i,:),testCont);
                            end
                            
                            ci95 = quantile(p_hat,[0.025 0.975],3);
                            cols = get(gca,'ColorOrder');
                            for r = 1:3
                                h=fill([diff(testCont,[],2);flipud(diff(testCont,[],2))],[ci95(:,r,1);flipud(ci95(:,r,2))],'b');
                                h.EdgeColor=[1 1 1];
                                h.FaceAlpha=0.5;
                                %                                 h.EdgeAlpha=h.FaceAlpha;
                                h.FaceColor=cols(r,:);
                                %                                 h.EdgeColor=cols(r,:);
                            end
                        catch
                            warning('Skipping plotting null laser effect bounds');
                        end
                        
                        if site == 1
                            ylabel(['ped: ' num2str(cVals(ped))]);
                        else
                            set(gca,'ytick','','ycolor','w');
                        end
                        
                        if ped ~= numPedestals
                            set(gca,'xtick','','xcolor','w');
                        end
                        
                        if ped == 1
                            title(siteLabels{site})
                        end
                        
                        %Overlay laser predictions
                        if length(p_L)==2
                            p_L = [p_L 0 0];
                        end
                        p_Total = p_nL + p_L;
                        p_hat = obj.calculatePhat(ZL_nL,ZR_nL,zeros(length(testCont),4),p_Total,testCont);
                        set(gca,'ColorOrderIndex',1);
                        plot(diff(testCont,[],2),p_hat,'linewidth',3);
                        xlim([-1 1]*(max(cVals)+0.1)); ylim([0 1]);
                    end
                    
                end
                
                set(gcf,'color','w');
            end
            
        end
        
        function plotData(obj,what,varargin)
            
            if cell2mat(strfind(varargin,'withSig')) == 1
                withSig = 1;
            else
                withSig = 0;
            end
            
            numSubjects = length(obj.names);
            val_stim_allSubj = cell(1,numSubjects);
            val_resp_allSubj = cell(1,numSubjects);
            for n = 1:numSubjects
                %Get metric of interest for different stimuli and responses
                figString = [obj.expt ' ' obj.names{n} ' delta ' what];
                figure('name',figString,'position',[300 500-50*n 1060 550]);
                
                %                 %corner conditions
                %                 maxC = max(obj.data{n}.contrast_cond);
                %                 c_0 = sum(obj.data{n}.contrast_cond,2)==0;
                %                 c_L = obj.data{n}.contrast_cond(:,1) == maxC(1) & obj.data{n}.contrast_cond(:,2) == 0;
                %                 c_R = obj.data{n}.contrast_cond(:,2) == maxC(2) & obj.data{n}.contrast_cond(:,1) == 0;
                %                 c_LR = obj.data{n}.contrast_cond(:,1) == maxC(1) & obj.data{n}.contrast_cond(:,2) == maxC(2);
                %                 stim = c_0 + 2*c_L + 3*c_R + 4*c_LR;
                %                 stim_labels = {'C=[0 0]','High CL','High CR','High CL&CR'};
                
                % % %                 %quadrant conditions
                %                 c_0 = sum(obj.data{n}.contrast_cond,2)==0;
                %                 c_L = obj.data{n}.contrast_cond(:,1) > obj.data{n}.contrast_cond(:,2);
                %                 c_R = obj.data{n}.contrast_cond(:,2) > obj.data{n}.contrast_cond(:,1);
                %                 c_LR = ( obj.data{n}.contrast_cond(:,1) == obj.data{n}.contrast_cond(:,2) ) & ~c_0;
                %                 stim = c_0 + 2*c_L + 3*c_R + 4*c_LR;
                %                 stim_labels = {'C=[0 0]','CL > CR','CL < CR','CL = CR'};
                % %
                %Detection + CL=CR conditions
                c_0 = sum(obj.data{n}.contrast_cond,2)==0;
                c_L = obj.data{n}.contrast_cond(:,1) > 0 & obj.data{n}.contrast_cond(:,2)==0;
                c_R = obj.data{n}.contrast_cond(:,2) > 0 & obj.data{n}.contrast_cond(:,1)==0;
                c_LR = ( obj.data{n}.contrast_cond(:,1) == obj.data{n}.contrast_cond(:,2) ) & ~c_0;
                stim = c_0 + 2*c_L + 3*c_R + 4*c_LR;
                stim_labels = {'C=[0 0]','CL only','CR only','CL = CR'};
                
                resp = obj.data{n}.response;
                resp_labels = {'Chose Left','Chose Right','Chose NoGo'};
                laser = obj.data{n}.laserIdx;
                
                h=uicontrol('Style','text','String', figString,'Units','normalized','Position', [0 0.9 1 0.1]);
                set(h,'Foregroundcolor','r','FontSize',20,'Fontname','Helvetica','Fontweight','bold');
                switch(what)
                    case 'performance'
                        val = (obj.data{n}.feedbackType==1);
                    case 'RT'
                        %                         numSessions = max(obj.data{n}.sessionID);
                        %                         rt_z = zeros(numSessions,1);
                        %                         for sessionID = 1:numSessions %Z-score the RTs for each session and each choice
                        %                             for choice = 1:2
                        %                                 idx = obj.data{n}.sessionID==sessionID & resp==choice;
                        %                                 rt_z(idx) = zscore(obj.data{n}.RT(idx));
                        %                             end
                        %                         end
                        %                         val = rt_z;
                        
                        val = obj.data{n}.RT;
                    case 'chooseL'
                        val = (resp==1);
                    case 'chooseR'
                        val = (resp==2);
                    case 'chooseNG'
                        val = (resp==3);
                    case 'pupil'
                        val = obj.data{n}.pupilenergy;
                    otherwise
                        %simulate data from model fits
                        [ZL,ZR] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                        offset = obj.fitData.nonLaserParams{n}(obj.data{n}.sessionID,:);
                        cont = obj.data{n}.contrast_cond;
                        resp = nan(length(cont),1);
                        for t = 1:size(cont,1)
                            if laser(t)==0
                                p = obj.calculatePhat(ZL,ZR,offset(t,:),[0 0 0 0],cont(t,:));
                            else
                                p = obj.calculatePhat(ZL,ZR,offset(t,:),obj.fitData.params{n}(laser(t),:),cont(t,:));
                            end
                            resp(t,1) = find(mnrnd(1,p));
                        end
                        switch(what)
                            case 'chooseL_simulated'
                                val = (resp==1);
                            case 'chooseNG_simulated'
                                val = (resp==3);
                            case 'chooseR_simulated'
                                val = (resp==2);
                        end
                end
                
                val(stim==0) = [];
                resp(stim==0) = [];
                laser(stim==0) = [];
                stim(stim==0)= [];
                
                %Compute mean value over all laser sites, for each stimulus
                %or response value separately
                val_stim = pivottable(laser,stim,val,'nanmean');
                val_resp = pivottable(laser,resp,val,'nanmean');
                
                %compute change in metric from nonLaser condition
                val_stim = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:));
                val_resp = bsxfun(@minus,val_resp(2:end,:),val_resp(1,:));
                
                if withSig == 1
                    %Non-laser sampling to get null distribution of metric
                    NUM_SHUFFLES = 500;
                    val_stim_null = nan(NUM_SHUFFLES,4);
                    val_resp_null = nan(NUM_SHUFFLES,3);
                    tab = tabulate(laser);
                    prop = mean(tab(2:end,3)/100); %get avg proportion of laser trials
                    nLidx = laser==0; %indices of all non-laser trials, to then sample a subset from
                    
                    %On each iteration, sample from the non-laser trials. Then
                    %do the exact same computation on the metric of interest
                    for shuff = 1:NUM_SHUFFLES
                        disp(shuff);
                        laserID = laser;
                        fakeSiteTrials = randsample(find(nLidx),round(prop*sum(nLidx)));
                        laserID(fakeSiteTrials) = max(laserID)+1;
                        
                        shuf_stim = pivottable(laserID,stim,val,'mean');
                        shuf_resp = pivottable(laserID,resp,val,'mean');
                        val_stim_null(shuff,:) = shuf_stim(end,:) - shuf_stim(1,:);
                        val_resp_null(shuff,:) = shuf_resp(end,:) - shuf_resp(1,:);
                    end
                    
                    %compute p value for metric at each
                    %site/response/choice config
                    val_resp_pvalue = nan(size(val_resp,1),3);
                    val_stim_pvalue = nan(size(val_stim,1),4);
                    for t = 1:size(val_resp,1)
                        lessThan = mean(bsxfun(@le,val_resp(t,:),val_resp_null),1);
                        moreThan = mean(bsxfun(@ge,val_resp(t,:),val_resp_null),1);
                        val_resp_pvalue(t,:) = min([lessThan; moreThan],[],1);
                        
                        lessThan = mean(bsxfun(@le,val_stim(t,:),val_stim_null),1);
                        moreThan = mean(bsxfun(@ge,val_stim(t,:),val_stim_null),1);
                        val_stim_pvalue(t,:) = min([lessThan; moreThan],[],1);
                    end
                    
                else
                    val_resp_pvalue = zeros(size(val_resp,1),3);
                    val_stim_pvalue = zeros(size(val_stim,1),4);
                end
                
                if size(val_stim,2)==3 %if detection task
                    val_stim = [val_stim, zeros(size(val_stim,1),1)];
                end
                
                ap = obj.inactivationCoords{n}(:,1);
                ml = obj.inactivationCoords{n}(:,2);
                if obj.bilateral_flag==1
                    ap = [ap; ap];
                    ml = [ml; -ml];
                    val_stim = [val_stim; val_stim];
                    val_resp = [val_resp; val_resp];
                    val_resp_pvalue = [val_resp_pvalue;val_resp_pvalue];
                    val_stim_pvalue = [val_stim_pvalue;val_stim_pvalue];
                end
                
                %Plot metric
                for s = 1:4
                    subplot(2,4,s); %Separated by stimulus
                    out = val_stim(:,s);
                    dotSize = 25 + (val_stim_pvalue(:,s)<0.05)*75 + (val_stim_pvalue(:,s)<0.01)*100;
                    scatter(ml,ap,dotSize,out,'s','filled'); axis equal;
                    xlabel(stim_labels{s});
                    %                     caxis([-1 1]*max(abs(val_stim(:,s))));
                    caxis([-1 1]*max(abs(quantile(val_stim(:),[0.025 0.975]))));
                    colorbar;
                end
                
                for r = 1:3
                    subplot(2,3,r+3); %Separated by response
                    out = val_resp(:,r);
                    dotSize = 25 + (val_resp_pvalue(:,r)<0.05)*140 + (val_resp_pvalue(:,r)<0.01)*165;
                    scatter(ml,ap,dotSize,out,'s','filled'); axis equal;
                    xlabel(resp_labels{r});
                    %                     caxis([-1 1]*max(abs(val_resp(:,r))));
                    caxis([-1 1]*max(abs(quantile(val_resp(:),[0.025 0.975]))));
                    colorbar;
                end
                
                cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                colormap(flipud(cmap));
                %
                set(findobj(get(gcf,'children'),'Type','Axes'),'xtick','','ytick','','box','off','xcolor','k','ycolor','w');
                set(gcf,'Color','w');
                
                val_stim_allSubj{n} = val_stim;
                val_resp_allSubj{n} = val_resp;
                drawnow;
            end
            
            
            %Plot average over mice
            ave_stim = mean(cat(3,val_stim_allSubj{:}),3);
            ave_resp = mean(cat(3,val_resp_allSubj{:}),3);
            figString = [obj.expt ' Average delta ' what];
            figure('name',figString,'position',[1500 200 1060 550]);
            h=uicontrol('Style','text','String', figString,'Units','normalized','Position', [0 0.9 1 0.1]);
            set(h,'Foregroundcolor','r','FontSize',20,'Fontname','Helvetica','Fontweight','bold');
            for s = 1:4
                subplot(2,4,s); %Separated by stimulus
                out = ave_stim(:,s);
                scatter(ml,ap,200,out,'s','filled'); axis equal;
                xlabel(stim_labels{s});
                %                 caxis([-1 1]*max(abs(ave_stim(:,s))));
                caxis([-1 1]*max(abs(quantile(ave_stim(:),[0.025 0.975]))));
                colorbar;
            end
            
            for r = 1:3
                subplot(2,3,r+3); %Separated by response
                out = ave_resp(:,r);
                scatter(ml,ap,380,out,'s','filled'); axis equal;
                xlabel(resp_labels{r});
                %                 caxis([-1 1]*max(abs(ave_resp(:,r))));
                caxis([-1 1]*max(abs(quantile(ave_resp(:),[0.025 0.975]))));
                colorbar;
            end
            
            cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
            colormap(flipud(cmap));
            set(findobj(get(gcf,'children'),'Type','Axes'),'xtick','','ytick','','box','off','xcolor','k','ycolor','w');
            set(gcf,'Color','w');
        end
        
        function plotData_RThist(obj)
            numSubjects = length(obj.names);
            choiceLabel = {'Left choices','Right choices'};
            for n = 1:numSubjects
                for r = 1:2
                    figString = [obj.names{n} ' ' choiceLabel{r}];
                    figure('name',figString,'color','w','WindowScrollWheelFcn',@callback_RTHIST);
                    h=uicontrol('Style','text','String', figString,'Units','normalized','Position', [0 0.9 1 0.1]);
                    set(h,'Foregroundcolor','r','FontSize',20,'Fontname','Helvetica','Fontweight','bold');
                
                    nonLaserRT = obj.data{n}.RT(obj.data{n}.laserIdx==0 & obj.data{n}.response==r);
                    for site = 1:size(obj.inactivationCoords{n},1);
                        
                        LaserRT = obj.data{n}.RT(obj.data{n}.laserIdx==site & obj.data{n}.response==r);
                        posIn = obj.inactivationCoords{n}(site,1:2)/10 + [0.58 0.465];
                        axes; hold on;
                        h1 = histogram(nonLaserRT,'BinMethod','scott','DisplayStyle','stairs','Normalization','pdf');
                        h2 = histogram(LaserRT,'BinMethod','scott','DisplayStyle','stairs','Normalization','pdf');
                        h1.LineWidth=1;
                        h2.LineWidth=1;
                        
                        set(gca,'Position',[posIn(2) posIn(1) 0.07 0.07],'box','off','ytick','');
                        set(gca,'xtick','','ycolor','w','xcolor','w');
                        xlim([0 1]);
                    end
                end
            end
            
            figure('name','POOLED TRIALS','color','w','WindowScrollWheelFcn',@callback_RTHIST);
            allRT = cellfun(@(s)(s.RT),obj.data,'uni',0)'; allRT = cat(1,allRT{:});
            allR = cellfun(@(s)(s.response),obj.data,'uni',0)'; allR = cat(1,allR{:});
            allL = cellfun(@(s)(s.laserIdx),obj.data,'uni',0)'; allL = cat(1,allL{:});
            nonLaserRT = allRT(allL==0 & allR==r);
            for site = 1:size(obj.inactivationCoords{1},1);
                
                LaserRT = allRT(allL==site & allR==r);
                posIn = obj.inactivationCoords{n}(site,1:2)/10 + [0.58 0.465];
                axes; hold on;
                h1 = histogram(nonLaserRT,'BinMethod','scott','DisplayStyle','stairs','Normalization','pdf');
                h2 = histogram(LaserRT,'BinMethod','scott','DisplayStyle','stairs','Normalization','pdf');
                h1.LineWidth=1;
                h2.LineWidth=1;
                
                set(gca,'Position',[posIn(2) posIn(1) 0.07 0.07],'box','off','ytick','');
                set(gca,'xtick','','ycolor','w','xcolor','w');
                xlim([0 1]);
            end
            
        end
        
        function plotFit_Comparison(obj)
        %Compare model-free delta pL,pR,pNG with model predicted delta pL,pR,pNG
            numSubjects = length(obj.names);
            for n = 1:numSubjects
                figure('name',obj.names{n},'color','w');
                cont = obj.data{n}.contrast_cond;
                c_0 = sum(cont,2)==0;
                c_L = cont(:,1) > 0 & cont(:,2)==0;
                c_R = cont(:,2) > 0 & cont(:,1)==0;
                c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
                stim = c_0 + 2*c_L + 3*c_R + 4*c_LR;
                stim_labels = {'C=[0 0]','CL only','CR only','CL = CR'};
                row_labels = {'dL','dR','dNG'};
                resp = obj.data{n}.response;
                laser = obj.data{n}.laserIdx;
                
                resp(stim==0) = [];
                laser(stim==0) = [];
                cont(stim==0,:) = [];
                stim(stim==0)= [];
                
                AX = [];
                %Plot actual data
                for r = 1:3
                    val_stim = pivottable(laser,stim,resp==r,'mean');
                    val_stim = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:));
                    
                    ap = obj.inactivationCoords{n}(:,1);
                    ml = obj.inactivationCoords{n}(:,2);
                    if obj.bilateral_flag==1
                        ap = [ap; ap];
                        ml = [ml; -ml];
                        val_stim = [val_stim; val_stim];
                    end
                    
                    for s = 1:4
                        AX(r,s)=subplot(3,4,s + 4*(r-1));
                        scatter(ml,ap,200,val_stim(:,s),'s','filled'); axis equal;
                        caxis([-1 1]*max(abs(val_stim(:))));
                        set(gca,'xtick','','ytick','','xcolor','w','ycolor','w');
                        hold on;
                        if r == 1
                            title('Actual & Pred');
                        end
                        
                        if s == 1
                            set(gca,'ycolor','k'); ylabel(row_labels{r});
                        end
                        
                        if r == 3
                            set(gca,'xcolor','k'); xlabel(stim_labels{s});
                        end
                    end
                end
                
                %Plot predictions
                nonlaserparams = mean(obj.fitData.nonLaserParams{n},1);
                [ZL,ZR,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                for s = 1:4
                    c = cont(stim==s,:);
                    nonLaserP = mean(obj.calculatePhat(ZL,ZR,nonlaserparams,[0 0 0 0],c),1);
                    LaserP = [];
                    for site = 1:size(obj.inactivationCoords{n},1)
                        LaserP(site,:) = mean(obj.calculatePhat(ZL,ZR,nonlaserparams + obj.fitData.params{n}(site,:),[0 0 0 0],c),1);
                    end
                    
                    dP = bsxfun(@minus,LaserP,nonLaserP);
                    
                    ap = obj.inactivationCoords{n}(:,1);
                    ml = obj.inactivationCoords{n}(:,2);
                    if obj.bilateral_flag==1
                        ap = [ap; ap];
                        ml = [ml; -ml];
                        dP = [dP; dP];
                    end
                    
                    for r = 1:3
                        scatter(AX(r,s),ml+11,ap,200,dP(:,r),'s','filled'); axis equal;
                    end
                end

                cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                colormap(flipud(cmap));
            end
        end
                
    end
    
    methods (Access=private)
        function out = cfn(obj,c,subj)
            out = (c.^obj.cfn_parameters{subj}(1))./(c.^obj.cfn_parameters{subj}(1) + obj.cfn_parameters{subj}(2).^obj.cfn_parameters{subj}(1));
        end
        
        function pupilenergy = addPupilData(obj,n,expt)
            
            try
                load(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' obj.names{n} '_PUPIL.mat']);
            catch
                pupilWindow = cell(1,length(obj.expRefs{n}));
                for b = 1:length(obj.expRefs{n})
                    try
                        eRef = obj.expRefs{n}{b};
                        vidFilename = ['\\zserver\Data\EyeCamera\' eRef(14:end) '\' eRef(1:10) '\' eRef(12) '\eye.mj2'];
                        v = VideoReader(vidFilename);
                        v.CurrentTime = 10;
                        imagesc(v.readFrame);
                        pupilWindow{b} = round(getrect);
                    catch
                    end
                end
                
                
                pupilenergy = cell(1,length(obj.expRefs{n}));
                for b = 1:length(obj.expRefs{n})
                    eRef = obj.expRefs{n}{b};
                    block = dat.loadBlock(eRef);
                    trials = block.trial;
                    
                    disp('Loading eye movie...');
                    pupilenergy{b} = nan(block.numCompletedTrials,1);
                    try
                        vidFilename = ['\\zserver\Data\EyeCamera\' eRef(14:end) '\' eRef(1:10) '\' eRef(12) '\eye.mj2'];
                        v = VideoReader(vidFilename);
                        
                        rect = pupilWindow{b};
                        rows = rect(2):(rect(2)+rect(4));
                        cols = rect(1):(rect(1)+rect(3));
                        
                        for t=1:block.numCompletedTrials
                            %                             disp(t);
                            start = trials(t).onsetToneSoundPlayedTime(1);
                            finish = trials(t).trialEndedTime;
                            v.CurrentTime = start;
                            a=0;
                            pixCounts=0;
                            while v.CurrentTime < finish
                                frame = 255 - v.readFrame;
                                frame = frame(rows,cols);
                                frame = frame>180;
                                imagesc(frame); title([num2str(t) '/' num2str(block.numCompletedTrials)]); drawnow;
                                pixCounts = pixCounts + sqrt(sum(frame(:)));
                                a=a+1;
                                v.CurrentTime = v.CurrentTime + 0.2;
                            end
                            
                            pupilenergy{b}(t,1) = pixCounts/a;
                        end
                        
                        %zscore pupil energy so that pupil energy is
                        %standardised between sessions
                        pupilenergy{b} = zscore(pupilenergy{b});
                    catch
                        warning('No eye data found.. Setting to NaNs');
                    end
                end
                save(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' obj.names{n} '_PUPIL.mat'],'pupilenergy');
            end
        end
        
        function pvalue = sig(obj)
            %use KDE to ask how likely an observed fitted point is under
            %the null
            numSubjects = length(obj.names);
            pvalue = cell(1,numSubjects);
            
            for n = 1:numSubjects
                p_null = obj.fitData.paramsNull{n};
                p = obj.fitData.params{n};
                
                %                 keyboard;
                pvalue{n} = nan(size(p));
                for i = 1:size(p,1)
                    for j = 1:size(p,2)
                        a = min([mean(p_null(:,j)>p(i,j)) mean(p_null(:,j)<p(i,j))]);
                        b = min([mean(p_null(:,j)>-p(i,j)) mean(p_null(:,j)<-p(i,j))]);
                        pvalue{n}(i,j) = a+b;
                    end
                end
                
                %TODO: specified parameter significance with the new
                %parameterisation
                
                %Bias params
                %                 ci_95_bL = quantile(p_null(:,1)./p_null(:,3),[0.025 0.975]);
                %                 ci_95_bR = quantile(p_null(:,2)./p_null(:,4),[0.025 0.975]);
                %                 sig_bL = ((p(:,1)./p(:,3)) < ci_95_bL(1)) | ((p(:,1)./p(:,3)) > ci_95_bL(2));
                %                 sig_bR = ((p(:,2)./p(:,4)) < ci_95_bR(1)) | ((p(:,2)./p(:,4)) > ci_95_bR(2));
                %
                %                 %Bias params
                %                 ci_95_bL = quantile(p_null(:,1),[0.025 0.975]);
                %                 ci_95_bR = quantile(p_null(:,2),[0.025 0.975]);
                %                 sig_bL = ((p(:,1)) < ci_95_bL(1)) | ((p(:,1)) > ci_95_bL(2));
                %                 sig_bR = ((p(:,2)) < ci_95_bR(1)) | ((p(:,2)) > ci_95_bR(2));
                %
                %                 %Sens params
                %                 ci_95_sL = quantile(p_null(:,3),[0.025 0.975]);
                %                 ci_95_sR = quantile(p_null(:,4),[0.025 0.975]);
                %                 sig_sL = (p(:,3) < ci_95_sL(1)) | (p(:,3) > ci_95_sL(2));
                %                 sig_sR = (p(:,4) < ci_95_sR(1)) | (p(:,4) > ci_95_sR(2));
                %
                %                 %Non-specific effect
                %                 s = mahal(p,p_null);
                %                 s = s/max(s);
                %                 sig_bL = s;
                %                 sig_bR = s;
                %                 sig_sL = s;
                %                 sig_sR = s;
                
                %non-specific effect test
                %Slice data along p-[0 0] vector and test against the
                %1D histogram
                %                 f = figure;
                %                 sigFlag = zeros(size(p,1),1);
                %                 for site = 1:size(p,1)
                %                     vec = p(site,:);
                %                     %perpendicular subspace
                %                     perp = null(vec);
                %                     %project null onto perp
                %                     proj = p_null*perp;
                %                     %isolate all null datapoints within a bound
                %                     distance = sqrt(sum(proj.^2,2));
                %                     idx = distance < 0.25;
                %                     %Project those isolated data onto the connecting vector
                %                     proj2 = p_null(idx,:)*vec';
                %
                %                     %is the projection of the parameter estimate on this
                %                     %vector outside the 95%ci of this null projection?
                %                     ci_95 = quantile(proj2,[0.025 0.975]);
                %                     ci_99 = quantile(proj2,[0.005 0.995]);
                %                     projParam = vec*vec';
                %
                %                     sigFlag(site) = min([mean(projParam>proj2) mean(projParam<proj2)]); %p value
                %                 end
                %
                %                 pvalue{n} = double(sigFlag);
            end
        end
        
        function X = getNonLaserDesignMatrix(obj,subj,transformed_contrasts,sessionIDs)
            numSessions = max(obj.data{subj}.sessionID);
            X = zeros(size(transformed_contrasts,1),3*numSessions);
            for t = 1:size(X,1)
                sid = sessionIDs(t);
                c = transformed_contrasts(t,:);
                X(t,3*sid - 2) = 1;
                X(t,3*sid - 1) = c(1);
                X(t,3*sid - 0) = c(2);
            end
            
            X=sparse(X);
            
            %             imagesc(X); drawnow;
        end
        
        function eRefs = getExpRefs(~,name,expt)
            try
                load(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' name '_EXPREFS.mat']);
            catch
                warning('Data not found, recalculating now...');
                
                eRefs = vertcat(dat.listExps(name));
                good = ones(length(eRefs),1);
                
                for b = 1:length(eRefs)
                    
                    try
                        p = load(dat.expFilePath(eRefs{b},'parameters','m'));
                        l = laserGLM(eRefs{b}); %check if any lasermanip problems
                        
                        switch(expt)
                            case 'sparse_bilateral_2D'
                                if p.parameters.numRepeats(23)~=1500
                                    good(b)=0;
                                end
                                
                                if p.parameters.rewardOnStimulus(2,end)~=2.5
                                    good(b)=0;
                                end
                                
                                block = dat.loadBlock(eRefs{b});
                                if block.numCompletedTrials<150
                                    good(b)=0;
                                end
                                
                                stim = p.parameters.visCueContrast;
                                if ~any(min(stim,[],1)>0)
                                    good(b)=0;
                                end
                                
                                if size(l.inactivationSite,1) ~= 26
                                    good(b)=0;
                                end
                                
                                
                            case 'sparse_unilateral_1D'
                                stim = p.parameters.visCueContrast;
                                if any(min(stim,[],1)>0)
                                    good(b)=0;
                                end
                                
                                %                 p.parameters.rewardOnStimulus
                                % max(p.parameters.rewardOnStimulus(2,:))
                                %                 if max(p.parameters.rewardOnStimulus(2,:)) ~= 1
                                %                     good(b)=0;
                                %                 end
                                %
                                l = laserGLM(eRefs{b});
                                if length(l.data.response)<150
                                    good(b)=0;
                                end
                                
                                if size(l.inactivationSite,1) < 5
                                    good(b)=0;
                                end
                                
                                
                                
                            case 'sparse_unilateral_2D'
                                stim = p.parameters.visCueContrast;
                                if ~any(min(stim,[],1)>0)
                                    good(b)=0;
                                end
                                
                                %                 p.parameters.rewardOnStimulus
                                % max(p.parameters.rewardOnStimulus(2,:))
                                %                 if max(p.parameters.rewardOnStimulus(2,:)) ~= 1
                                %                     good(b)=0;
                                %                 end
                                %
                                if length(l.data.response)<150
                                    good(b)=0;
                                end
                                
                                if size(l.inactivationSite,1) < 45
                                    good(b)=0;
                                end
                                
                                if ~any(l.data.laser(:,2)<0)
                                    good(b)=0;
                                end
                                
                            case 'sparse_unilateral_1D_highRes'
                                stim = p.parameters.visCueContrast;
                                if any(min(stim,[],1)>0)
                                    good(b)=0;
                                end
                                
                                if length(l.data.response)<150
                                    good(b)=0;
                                end
                                
                                if size(l.inactivationSite,1) < 5
                                    good(b)=0;
                                end
                                
                                if size(l.inactivationSite,1) > 10
                                    good(b)=0;
                                end
                                
                            case 'sparse_unilateral_2D_highRes'
                                
                                good(b)=0;
                                if strcmp(l.sessionLabel,'LineAlongLeftPPC')
                                    good(b)=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightPPC')
                                    good(b)=1;
                                elseif strcmp(l.sessionLabel,'LineAlongLeftVisual')
                                    good(b)=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightVisual')
                                    good(b)=1;
                                elseif strcmp(l.sessionLabel,'LineAlongLeftM2')
                                    good(b)=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightM2')
                                    good(b)=1;
                                end
                                
                            case 'sparse_unilateral_2D_PPC'
                                
                                good(b)=0;
                                if strcmp(l.sessionLabel,'LineAlongLeftPPC')
                                    good(b)=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightPPC')
                                    good(b)=1;
                                end
                                
                            case 'sparse_unilateral_2D_VIS'
                                good(b)=0;
                                
                                if strcmp(l.sessionLabel,'LineAlongLeftVisual')
                                    good(b)=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightVisual')
                                    good(b)=1;
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
                
                save(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' name '_EXPREFS.mat'],'eRefs');
            end
        end
        
        function logLik = calculateLogLik(obj,testParams,ZL,ZR,offset,inputs,responses)
            phat = obj.calculatePhat(ZL,ZR,offset,testParams,inputs);
            p = (responses==1).*phat(:,1) + (responses==2).*phat(:,2) + (responses==3).*phat(:,3);
            p(p<0)=0; %Correct for weird numerical error at extreme ends of probabilities;
            %             if any(p<0)
            %                 keyboard;
            %             end
            
            logLik = mean(log2(p));
        end
        
        function phat = calculatePhat(~,ZL,ZR,offset,testParams,inputs)
            zl = ZL(offset,testParams,inputs);
            zr = ZR(offset,testParams,inputs);
            pL = exp(zl)./(1+exp(zl)+exp(zr));
            pR = exp(zr)./(1+exp(zl)+exp(zr));
            pNG = 1 - pL - pR;
            
            phat = [pL pR pNG];
        end
        
    end
end