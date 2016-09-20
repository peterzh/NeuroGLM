classdef omnibusLaserGLM
    properties
        names = {'Spemann','Murphy','Morgan','Whipple'};
        expRefs;
        inactivationCoords;
        data;
        fitData;
        crossval;
        fitMethod='fmincon';
        guess_bpt = [];
    end
    
    properties (Access=public)
        cfn_parameters;
        moreData = 0;
        significance_flag = 0;
    end
    
    methods
        function obj = omnibusLaserGLM(expt)
            disp('Collating expRefs...');
            numSubjects = length(obj.names);
            obj.expRefs = cell(1,numSubjects);
            for n = 1:numSubjects
                obj.expRefs{n} = obj.getExpRefs(obj.names{n},expt);
            end
            
            disp('Loading data...');
            obj.inactivationCoords = cell(1,numSubjects);
            obj.data = cell(1,numSubjects);
            obj.cfn_parameters = cell(1,numSubjects);
            for n = 1:numSubjects
                D = struct;
                c50_parameters = [];
                for session = 1:length(obj.expRefs{n})
                    d = laserGLM(obj.expRefs{n}{session});
                    
                    if obj.moreData==1
                        d = d.addData('lick');
                        d = d.addData('pupil');
                    end
                    d = d.data;
                    d = structfun(@(x)(x(6:(end-14),:)),d,'uni',0); %trim first 5 trials and last 15
                    
                    %Fit simple C50 models to estimate c50 & n params
                    e = getrow(d,d.laserIdx==0);
                    g=GLM(e).setModel('C50-subset').fit;
                    c50_parameters(session,:) = g.parameterFits(5:6);
                    
                    d.sessionID = ones(length(d.response),1)*session;
                    D = addstruct(D,d);
                end
                D = rmfield(D,'laserIdx');
                D = getrow(D,D.repeatNum==1);
                l = laserGLM(D);
                
                obj.inactivationCoords{n} = l.inactivationSite;
                obj.data{n} = l.data;
                
                tab = tabulate(D.response);
                tab = tab(:,3)/100;
                obj.guess_bpt(n)=sum(tab.*log2(tab));
                
                obj.cfn_parameters{n} = mean(c50_parameters,1);
            end
        end
        
        function [ZL,ZR,numP] = getModel(obj,model,subj)
            
            %             cfn = @(c)(c.^obj.cfn_parameters{subj}(1))./(c.^obj.cfn_parameters{subj}(1) + obj.cfn_parameters{subj}(2).^obj.cfn_parameters{subj}(1));
            cfn = @(c)(obj.cfn(c,subj));
            switch(model)
                case 'bias'
                    ZL = @(offL,p,c)(offL + p(1) );
                    ZR = @(offR,p,c)(offR + p(2) );
                    numP = 2;
                case 'sub_sens'
                    ZL = @(offL,p,c)(offL + p(1).*cfn(c(:,1)) );
                    ZR = @(offR,p,c)(offR + p(2).*cfn(c(:,2)) );
                    numP = 2;
                case 'bias+sub_sens'
                    ZL = @(offL,p,c)(offL + p(1) + p(3).*cfn(c(:,1)) );
                    ZR = @(offR,p,c)(offR + p(2) + p(4).*cfn(c(:,2)) );
                    numP = 4;
                case 'bias+sub_sens_cent'
                    c_mean = mean(cfn(obj.data{subj}.contrast_cond));
                    ZL = @(offL,p,c)(offL + p(1) + p(3).*(cfn(c(:,1))-c_mean(1)) );
                    ZR = @(offR,p,c)(offR + p(2) + p(4).*(cfn(c(:,2))-c_mean(2)) );
                    numP = 4;
                case 'bias+full_sens'
                    ZL = @(offL,p,c)(offL + p(1) + p(3).*cfn(c(:,1)) + p(5).*cfn(c(:,2)) );
                    ZR = @(offR,p,c)(offR + p(2) + p(4).*cfn(c(:,2)) + p(6).*cfn(c(:,1)));
                    numP = 6;
            end
        end
        
        function obj = fit(obj,lambda)
            f=figure;
            obj.fitData.nonLaserModel = 'bias+sub_sens';
            obj.fitData.laserModel = 'bias';
            obj.fitData.lambda = lambda;
            
            numSubjects = length(obj.names);
            %             numSubjects = 1;
            obj.fitData.params = cell(1,numSubjects);
            obj.fitData.nonLaserParams = cell(1,numSubjects);
            for n = 1:numSubjects
                
                %First fit non-laser portion of the data
                
                %Manual optim method per-session
                [ZL_nL,ZR_nL,numP] = obj.getModel(obj.fitData.nonLaserModel,n);
                nL = obj.data{n}.laserIdx==0;
                
                sessions = unique(obj.data{n}.sessionID(nL));
                p_nL = nan(length(sessions),numP);
                for s = 1:length(sessions)
                    cont = obj.data{n}.contrast_cond(nL & obj.data{n}.sessionID==sessions(s),:);
                    resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                    
                    offset = zeros(length(resp),2);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS.^2) );
                    
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
                
%                 keyboard;
                %
                %                 %GLMNET METHOD:
                %                 nL = obj.data{n}.laserIdx==0;
                %                 cont = obj.data{n}.contrast_cond(nL,:);
                %                 resp = obj.data{n}.response(nL);
                %                 X = obj.getNonLaserDesignMatrix(n,obj.cfn(cont,n),obj.data{n}.sessionID(nL));
                %                 opts = []; opts.intr=0; %exclude bias
                %                 fit = cvglmnet(X,resp,'multinomial',glmnetSet(opts));
                
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
                    offset = [];
                    for t = 1:sum(L_idx)
                        offset(t,:) = [ZL_nL(0,p_nL(sess(t),:),cont(t,:)),...
                                       ZR_nL(0,p_nL(sess(t),:),cont(t,:)) ];
                    end
                    
                    %Offset from glmnet method
                    %                     X = obj.getNonLaserDesignMatrix(n,obj.cfn(cont,n),obj.data{n}.sessionID(L_idx));
                    %                     offset = cvglmnetPredict(fit,X,obj.fitData.lambda,'link'); offset=[offset(:,1)-offset(:,3) offset(:,2)-offset(:,3)];
                    
                    %Use that offset to manually fit a model at each
                    %inactivation site
                    [ZL,ZR,numP] = obj.getModel(obj.fitData.laserModel,n);
                    objective = @(PARAMETERS) ( -obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS.^2) );
                    %                     keyboard;
%                     LB = ones(1,numP)*-10;
%                     UB = ones(1,numP)*10;
%                     
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
                obj.fitData.paramsLL{n} = -loglik;
            end
             
            disp('Done!'); close(f);
        end
        
        function obj = withBootstrap(obj)
            %Fits the non-laser model under different levels of resampling
            %Just use GLMNET fitting as that's quicker to code
            numSubjects = length(obj.names);
%             figure;
            obj.fitData.paramsSTD = cell(1,numSubjects);
            for n = 1:numSubjects
                [ZL_nL,ZR_nL,numP] = obj.getModel(obj.fitData.nonLaserModel,n);
                
                numSessions = max(obj.data{n}.sessionID);
                for s = 1:numSessions
                    nL = obj.data{n}.laserIdx==0 & obj.data{n}.sessionID==s;
                    
                    p_nL = [];
                    for iter = 1:100
                        nL_bootstrap = randsample(find(nL),sum(nL),true);
                        cont = obj.data{n}.contrast_cond(nL_bootstrap,:);
                        resp = obj.data{n}.response(nL_bootstrap);
                        
                        offset = zeros(length(resp),2);
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
                    
                    obj.fitData.nonLaserParamsSTD{n}(s,:) = std(p_nL);
%                     boxplot(p_nL{n});
                end
            end
        end
%         
        function obj = fitNull(obj)
            figure;
            kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7); hold on;
            cols = [linspace(0,1,52)' linspace(1,0,52)' zeros(52,1) ];
            cols(obj.inactivationCoords{1}(:,2)>0,3) = 1;
            scatter(obj.inactivationCoords{1}(:,2),obj.inactivationCoords{1}(:,1),200,cols,'filled');
            
            
            numSubjects = length(obj.names);
            %             numSubjects=1;
            numIter = 1000;
            
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
%                     plot(obj.fitData.params{n}(:,1)./obj.fitData.params{n}(:,3),obj.fitData.params{n}(:,3),'r.','markersize',15);
                elseif numLP == 2
%                     plot(obj.fitData.params{n}(:,1),obj.fitData.params{n}(:,2),'r.','markersize',15);
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
                    [ZL_nL,ZR_nL,numP] = obj.getModel(obj.fitData.nonLaserModel,n);
                    
                    sessions = unique(obj.data{n}.sessionID(nL));
                    p_nL = nan(length(sessions),numP);
                    for s = 1:length(sessions)
                        cont = obj.data{n}.contrast_cond(nL & obj.data{n}.sessionID==sessions(s),:);
                        resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                        
                        offset = zeros(length(resp),2);
                        objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS.^2));
                        
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
                    
                    %GLMNET METHOD
                    %                     XnL = X; XnL(q_idx,:) = [];
                    %                     RnL = resp; RnL(q_idx) = [];
                    %                     opts = []; opts.intr=1; %include global bias
                    %                     fit = cvglmnet(XnL,RnL,'multinomial',glmnetSet(opts));
                    
                    %Now fit the laser site using the offset defined from
                    %the non-laser component
%                     XL = X(q_idx,:);
%                     RL = resp(q_idx);
%                     contL = cont(q_idx,:);
                    
%                     offset = cvglmnetPredict(fit,XL,lambda,'link');
%                     offset=[offset(:,1)-offset(:,3) offset(:,2)-offset(:,3)];

                    cont = obj.data{n}.contrast_cond(q_idx,:);
                    resp = obj.data{n}.response(q_idx);
                    sess = obj.data{n}.sessionID(q_idx);
                    
                    %offset 
                    offset = [];
                    for t = 1:length(q_idx)
                        offset(t,:) = [ZL_nL(0,p_nL(sess(t),:),cont(t,:)),...
                                       ZR_nL(0,p_nL(sess(t),:),cont(t,:)) ];
                    end

                    [ZL,ZR,numP] = obj.getModel(obj.fitData.laserModel,n);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) + obj.fitData.lambda*sum(PARAMETERS.^2));
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                            [p,loglik(iter,1)] = Opt.solve;
                            obj.fitData.paramsNull{n}(iter,:) = p';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                            [p,loglik(iter,1)] = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                            obj.fitData.paramsNull{n}(iter,:) = p;
                    end
                    if numLP == 4
                        plot(p(1)/p(3),p(3),'k.','markersize',10); drawnow; xlabel('B_L/S_L'); ylabel('S_L');
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
        
        function plotNull(obj)
            figure;
            kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7); hold on;
            cols = [linspace(0,1,52)' linspace(1,0,52)' zeros(52,1) ];
            cols(obj.inactivationCoords{1}(:,2)>0,3) = 1;
%             cols(:,3) = obj.inactivationCoords{1}(:,2)/7 + .5;
            scatter(obj.inactivationCoords{1}(:,2),obj.inactivationCoords{1}(:,1),200,cols,'filled');
            
            numLP = size(obj.fitData.paramsNull{1},2);
            
            figure('name','Sampling from null distribution of laser parameters');
            numSubjects = length(obj.names);
            for n = 1:numSubjects
                
                if numLP == 4
                    
                    %Plot original B vs S
                    subplot(4,numSubjects,n); title(obj.names{n}); hold on;
                    plot(obj.fitData.paramsNull{n}(:,1),obj.fitData.paramsNull{n}(:,3),'.','markerfacecolor',[1 1 1]*0.5,'markeredgecolor',[1 1 1]*0.5);
                    scatter(obj.fitData.params{n}(:,1),obj.fitData.params{n}(:,3),30,cols,'filled');
                    xlabel('\Delta B_L'); ylabel('\Delta S_L'); axis square;
                    
                    subplot(4,numSubjects,numSubjects+n); hold on;
                    plot(obj.fitData.paramsNull{n}(:,2),obj.fitData.paramsNull{n}(:,4),'.','markerfacecolor',[1 1 1]*0.5,'markeredgecolor',[1 1 1]*0.5);
                    scatter(obj.fitData.params{n}(:,2),obj.fitData.params{n}(:,4),30,cols,'filled');
                    xlabel('\Delta B_R'); ylabel('\Delta S_R'); axis square;
                    
                    %Plot transformed bias effect for laser
%                     keyboard;
                    numTrials = length(obj.data{n}.sessionID);
                    P = nan(numTrials,4);
                    for t = 1:numTrials
                        sess = obj.data{n}.sessionID(t);                 
                        P(t,:) = obj.fitData.nonLaserParams{n}(sess,:);
                    end
                    ave_nonLaserP = mean(P);
                    
                    %add non-laser component to null, and to site
                    %parameters, in order to then transform to the new bias
                    %form
                    pnull_incNL = bsxfun(@plus,obj.fitData.paramsNull{n},ave_nonLaserP);
                    p_incNL = bsxfun(@plus,obj.fitData.params{n},ave_nonLaserP);
                    
                    subplot(4,numSubjects,2*numSubjects + n); hold on;
                    plot(pnull_incNL(:,1)./pnull_incNL(:,3),pnull_incNL(:,3),'.','markerfacecolor',[1 1 1]*0.5,'markeredgecolor',[1 1 1]*0.5);
                    scatter(p_incNL(:,1)./p_incNL(:,3),p_incNL(:,3),30,cols,'filled');
                    xlabel('B_{LnL} + B_{LL}/S_{LnL} + S_{LL}'); ylabel('S_{LnL} + S_{LL}'); axis square;
                    
                    subplot(4,numSubjects,3*numSubjects + n); hold on;
                    plot(pnull_incNL(:,2)./pnull_incNL(:,4),pnull_incNL(:,4),'.','markerfacecolor',[1 1 1]*0.5,'markeredgecolor',[1 1 1]*0.5);
                    scatter(p_incNL(:,2)./p_incNL(:,4),p_incNL(:,4),30,cols,'filled');
                    xlabel('B_{RnL} + B_{RL}/S_{RnL} + S_{RL}'); ylabel('S_{RnL} + S_{RL}'); axis square;
                    
                elseif numLP == 2
                    subplot(1,numSubjects,n); title(obj.names{n}); hold on;
                    plot(obj.fitData.paramsNull{n}(:,1),obj.fitData.paramsNull{n}(:,2),'.','markerfacecolor',[1 1 1]*0.5,'markeredgecolor',[1 1 1]*0.5);
                    scatter(obj.fitData.params{n}(:,1),obj.fitData.params{n}(:,2),30,cols,'filled');
                    xlabel('B_L'); ylabel('B_R'); axis square;
                end
            end
            
            set(gcf,'color','w');
        end
        
        function paramsSIG = sig(obj)
            %use KDE to ask how likely an observed fitted point is under
            %the null
            numSubjects = length(obj.names);
            paramsSIG = cell(1,numSubjects);
           
            for n = 1:numSubjects
                p_null = obj.fitData.paramsNull{n};
                p = obj.fitData.params{n};
                
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
                %Non-specific effect
                s = mahal(p,p_null);
                s = s/max(s);
                sig_bL = s;
                sig_bR = s;
                sig_sL = s;
                sig_sR = s;
                
                %Slice data along p-[0 0] vector and test against the
                %1D histogram
                for site = 1:size(p,1)
                    
                end
                
                paramsSIG{n} = double([sig_bL sig_bR sig_sL sig_sR]);
                
                
            end
        end
        
        function obj = fitCV(obj,numFolds)
            f=figure;
            obj.crossval = cell(1,length(obj.names));
            laserModels = {'sub_+bias','sub_+sens','sub_+bias+sens'};
            lambda = 0.001;
            %             lambda = 'lambda_min';
            numSubjects = length(obj.names);
            for n = 1:numSubjects
                
                cv = cvpartition(obj.data{n}.response,'KFold',numFolds);
                phat = nan(length(obj.data{n}.response),length(laserModels)+1);
                for fold = 1:cv.NumTestSets
                    
                    trainC = obj.data{n}.contrast_cond(cv.training(fold),:);
                    trainR = obj.data{n}.response(cv.training(fold));
                    trainL = obj.data{n}.laserIdx(cv.training(fold));
                    trainSID = obj.data{n}.sessionID(cv.training(fold));
                    
                    %Fit non-laser model to training set
                    nL_idx = trainL==0;
                    cont = trainC(nL_idx,:);
                    resp = trainR(nL_idx);
                    X = obj.getNonLaserDesignMatrix(n,obj.cfn(cont,n),trainSID(nL_idx));
                    opts = []; opts.intr=0; %exclude bias
                    fit = cvglmnet(X,resp,'multinomial',glmnetSet(opts));
                    
                    p_L = cell(size(obj.inactivationCoords{n},1),length(laserModels));
                    for site = 1:size(obj.inactivationCoords{n},1)
                        L_idx = trainL==site;
                        cont = trainC(L_idx,:);
                        resp = trainR(L_idx);
                        X = obj.getNonLaserDesignMatrix(n,obj.cfn(cont,n),trainSID(L_idx));
                        offset = cvglmnetPredict(fit,X,lambda,'link'); offset=[offset(:,1)-offset(:,3) offset(:,2)-offset(:,3)];
                        
                        for model = 1:length(laserModels)
                            [ZL,ZR,numP] = obj.getModel(laserModels{model},n);
                            objective = @(PARAMETERS) (obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp));
                            
                            switch(obj.fitMethod)
                                case 'opti'
                                    options = optiset('display','final','solver','NOMAD');
                                    Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                                    p_L{site,model} = Opt.solve';
                                    
                                case 'fmincon'
                                    options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                                    p_L{site,model} = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                            end
                        end
                    end
                    
                    %Now test these fits at the sites in the test set
                    testC = obj.data{n}.contrast_cond(cv.test(fold),:);
                    testR = obj.data{n}.response(cv.test(fold));
                    testL = obj.data{n}.laserIdx(cv.test(fold));
                    testSID = obj.data{n}.sessionID(cv.test(fold));
                    
                    for site = 1:size(obj.inactivationCoords{n},1)
                        L_idx = testL==site;
                        cont = testC(L_idx,:);
                        resp = testR(L_idx);
                        
                        if ~isempty(resp)
                            X = obj.getNonLaserDesignMatrix(n,obj.cfn(cont,n),testSID(L_idx));
                            offset = cvglmnetPredict(fit,X,lambda,'link'); offset=[offset(:,1)-offset(:,3) offset(:,2)-offset(:,3)];
                            
                            %non-laser model
                            ZL = @(offL,p,c)(offL);
                            ZR = @(offR,p,c)(offR);
                            p = obj.calculatePhat(ZL,ZR,offset,[],[]);
                            phat(cv.test(fold) & obj.data{n}.laserIdx==site,1) = (resp==1).*p(:,1) + (resp==2).*p(:,2) + (resp==3).*p(:,3);
                            
                            %now laser models
                            for model = 1:length(laserModels)
                                [ZL,ZR,~] = obj.getModel(laserModels{model},n);
                                p = obj.calculatePhat(ZL,ZR,offset,p_L{site,model},cont);
                                phat(cv.test(fold) & obj.data{n}.laserIdx==site,model+1) = (resp==1).*p(:,1) + (resp==2).*p(:,2) + (resp==3).*p(:,3);
                            end
                        end
                    end
                end
                
                for site = 1:size(obj.inactivationCoords{n},1)
                    for model = 1:size(phat,2)
                        obj.crossval{n}(site,model) = mean(log2(phat(obj.data{n}.laserIdx==site,model))) - obj.guess_bpt(n);
                    end
                end
            end
            
            disp('Done!');close(f);
        end
        
        function plotCrossval(obj)
            %             alpha=0.05; beta=0.6;
            %
            numSubjects = length(obj.names);
            model_labels = {'no laser terms','+bias','+sens','+bias+sens'};
            
            figure;
            for n = 1:numSubjects
                cols = [];
                
                %                 for m = 1:length(model_labels)
                %                     subplot(numSubjects,4,4*n-4 + m);
                %                     ll = obj.crossval{n};
                %                     ll(:,m)=[];
                %                     scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),300,obj.crossval{n}(:,m)-max(ll,[],2),'s','filled');
                %                     axis equal; caxis([-1 1]*0.1);
                %                     cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                %                         linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                %                     colormap(cmap);
                %                 end
                
                
                for site = 1:size(obj.inactivationCoords{n},1)
                    [y,winner] = max(obj.crossval{n}(site,:));
                    srt=sort(obj.crossval{n}(site,:));
                    
                    if (obj.crossval{n}(site,1)-y) > -0.05
                        %                         keyboard;
                        cols(site,:) = [1 1 1];
                    end
                    
                    
                    if winner == 2
                        cols(site,:) = [1 0 0];
                    elseif winner == 3
                        cols(site,:) = [0 1 0];
                    elseif winner == 4
                        cols(site,:) = [1 1 0];
                    end
                    
                    %                      cols(site,:) =  cols(site,:)*30*(y-srt(end-1));
                end
                subplot(1,numSubjects,n);
                scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),300,cols,'s','filled'); axis equal;
                title(obj.names{n});
            end
        end
        
        function obj = plotFit(obj)
            numSubjects = length(obj.names);
            
            %Plot session-by-session bias and sensitivity parameters
            figure('name','Session by session non-laser bias and sensitivity changes');
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
            set(gcf,'color','w');
            
            %Plot non-laser predictions against actual psychometric data.
%             figure('name','non-laser predictions against actual psychometric data');
            cols = {[0 0.4470 0.7410],...
                    [0.8500 0.3250 0.0980],...
                    [0.9290 0.6940 0.1250]};

            for n = 1:numSubjects
                figure('name',[obj.names{n} ' non-laser predictions against actual psychometric data']);
                [ZL_nL,ZR_nL,numP] = obj.getModel(obj.fitData.nonLaserModel,n);
                
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
                        p_hat = obj.calculatePhat(ZL_nL,ZR_nL,zeros(length(testCont),2),p_nL,testCont);
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
            
            %Plot laser parameter maps
            figure('name','Laser parameter maps');
            for n = 1:numSubjects
                p_L = obj.fitData.params{n};
                
                switch(obj.fitData.laserModel)
                    case 'bias'
                        plotVal = {p_L(:,1),p_L(:,2),p_L(:,1)-p_L(:,2), p_L(:,1)+p_L(:,2)};
                        param_labels = {'b_L','b_R','b_L-b_R','b_L+b_R'};
                    case 'bias+sub_sens'
                        plotVal = {p_L(:,1),p_L(:,2),p_L(:,3),p_L(:,4)};
                        param_labels = {'b_L','b_R','s_L','s_R'};

                end
                
                for i = 1:length(plotVal)
                    subplot(numSubjects,length(plotVal),length(plotVal)*n-length(plotVal) + i);
                    
                    if obj.significance_flag == 1
                        sigFlag = obj.sig;
                        
                        dotSize = 50 + sigFlag{n}(:,i)*200;
                        scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),dotSize,plotVal{i},'s','filled'); axis equal;

                    else
                        scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),200,plotVal{i},'s','filled'); axis equal;
                    end

                    caxis([-1 1]*5);
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
                colormap(cmap);
            end
            set(gcf,'color','w');

        end
        
        function modelFreePlot(obj,what)
            numSubjects = length(obj.names);
            bilateral_mirror = 0;
%             switch(what)
%                 case 'performance'
%                     figure('name','performance');
%                     kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
%                     clear PD;
%                     for s = 1:numSubjects
%                         subplot(1,numSubjects,s);
%                         imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
%                         set(imX,'alphadata',0.7);
%                         hold on;
%                         
%                         PF = 100*pivottable(obj.data{s}.laserIdx,[],obj.data{s}.feedbackType==1,'mean');
%                         
%                         % NULL TEST
%                         tab = tabulate(obj.data{s}.laserIdx);
%                         proportion = mean(tab(2:end,3)/100);
%                         NUM_SHUFFLES = 1000;
%                         nullP = nan(NUM_SHUFFLES,1);
%                         nonL_idx = obj.data{s}.laserIdx==0;
%                         
%                         for shuff = 1:NUM_SHUFFLES
%                             disp(shuff);
%                             L_idx = obj.data{s}.laserIdx;
%                             nL = randsample(find(nonL_idx),round(proportion*sum(nonL_idx)));
%                             L_idx(nL) = max(obj.data{s}.laserIdx)+1;
%                             
%                             shuf_p = 100*pivottable(L_idx,[],obj.data{s}.feedbackType==1,'mean');
%                             nullP(shuff) = shuf_p(end);
%                         end
%                         
%                         nullP_ci95 = quantile(nullP,[0.025 0.975]);
%                         nullP_ci99 = quantile(nullP,[0.005 0.995]);
%                         
%                         
%                         ap = obj.inactivationCoords{s}(:,1);
%                         ml = obj.inactivationCoords{s}(:,2);
%                         out = PF(2:end);
%                         if bilateral_mirror == 1
%                             ap = [ap;ap];
%                             ml = [ml;-ml];
%                             out = [out;out];
%                         end
%                         
%                         %sig testing
%                         sig95 = double((out < nullP_ci95(1)) | (nullP_ci95(2) < out));
%                         sig99 = double((out < nullP_ci99(1)) | (nullP_ci99(2) < out));
%                         sig = sum([sig95 sig99],2);
%                         sig(sig==0)=0.2;
%                         
%                         s1=scatter(ml,ap,sig*80,out,'s','filled'); axis equal; colorbar;
%                         %         hold on;
%                         s1.MarkerEdgeColor=[1 1 1]*0.9;
%                         
%                         %     s1=scatter(ml,ap,150,out,'s','filled'); axis equal; colorbar;
%                         ylim([-6 4]); %s1.MarkerEdgeColor=[1 1 1]*0.8;
%                         cax = PF(1) + [-1 1]*max(abs(PF(2:end)-PF(1)));
%                         try
%                             PD(:,s)=PF(2:end)-PF(1);
%                         catch
%                         end
%                         %     cax(cax>1)=1;
%                         caxis(cax);
%                         hold off;
%                         title(obj.names{s})
%                         xlim([-5 5]);
%                         
%                         cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
%                             linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
%                         colormap(flipud(cmap));
%                         
%                         set(gca,'xtick','','ytick','');
%                         set(gca,'box','off','XColor','w','YColor','w');
%                         
%                     end
%                     
%                     cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
%                         linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
%                     colormap(flipud(cmap));
%                     
%                     set(gca,'xtick','','ytick','');
%                     set(gca,'box','off','XColor','w','YColor','w');
%                     set(gcf,'Color','w');
%                     
%                 case 'RT'
%                 case 'choices'
%             end
%             
%             
            
            
            numSubjects = length(obj.names);
            val_stim_allSubj = cell(1,numSubjects);
            val_resp_allSubj = cell(1,numSubjects);
            for n = 1:numSubjects
                %Get metric of interest for different stimuli 
                figure('name',[obj.names{n} ' \Delta' what]);
                
%                 %corner conditions
%                 maxC = max(obj.data{n}.contrast_cond);
%                 c_0 = sum(obj.data{n}.contrast_cond,2)==0;
%                 c_L = obj.data{n}.contrast_cond(:,1) == maxC(1) & obj.data{n}.contrast_cond(:,2) == 0;
%                 c_R = obj.data{n}.contrast_cond(:,2) == maxC(2) & obj.data{n}.contrast_cond(:,1) == 0;
%                 c_LR = obj.data{n}.contrast_cond(:,1) == maxC(1) & obj.data{n}.contrast_cond(:,2) == maxC(2);
%                 stim = c_0 + 2*c_L + 3*c_R + 4*c_LR;
%                 stim_labels = {'C=[0 0]','High CL','High CR','High CL&CR'};
                
%                 %quadrant conditions
                c_0 = sum(obj.data{n}.contrast_cond,2)==0;
                c_L = obj.data{n}.contrast_cond(:,1) > obj.data{n}.contrast_cond(:,2);
                c_R = obj.data{n}.contrast_cond(:,2) > obj.data{n}.contrast_cond(:,1);
                c_LR = ( obj.data{n}.contrast_cond(:,1) == obj.data{n}.contrast_cond(:,2) ) & ~c_0;
                stim = c_0 + 2*c_L + 3*c_R + 4*c_LR;
                stim_labels = {'C=[0 0]','CL > CR','CL < CR','CL = CR'};
                
                resp = obj.data{n}.response;
                resp_labels = {'Chose Left','Chose Right','Chose NoGo'};
                                
                switch(what)
                    case 'performance'
                        val_stim = 100*pivottable(obj.data{n}.laserIdx,stim,obj.data{n}.feedbackType==1,'mean');
                        val_resp = 100*pivottable(obj.data{n}.laserIdx,resp,obj.data{n}.feedbackType==1,'mean');
                    case 'RT'
                        numSessions = max(obj.data{n}.sessionID);
                        rt_z = nan(numSessions,1);
                        for sessionID = 1:numSessions
                            for choice = 1:3
                                idx = obj.data{n}.sessionID==sessionID & resp==choice;
                                rt_z(idx) = zscore(obj.data{n}.RT(idx));
                            end
                        end
                        
                        val_resp = pivottable(obj.data{n}.laserIdx,resp,rt_z,'mean');
                        val_stim = pivottable(obj.data{n}.laserIdx,stim,rt_z,'mean');
                    case 'chooseL'
                        val_stim = 100*pivottable(obj.data{n}.laserIdx,stim,resp==1,'mean');
                        val_resp = 100*pivottable(obj.data{n}.laserIdx,resp,resp==1,'mean');
                    case 'chooseR'
                        val_stim = 100*pivottable(obj.data{n}.laserIdx,stim,resp==2,'mean');
                        val_resp = 100*pivottable(obj.data{n}.laserIdx,resp,resp==2,'mean');
                end
                
                %compute change in metric from nonLaser condition
                val_stim = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:)); 
                val_resp = bsxfun(@minus,val_resp(2:end,:),val_resp(1,:)); 
                
                val_stim_allSubj{n} = val_stim;
                val_resp_allSubj{n} = val_resp;
                
                %Plot metric
                for s = 1:4
                    subplot(2,4,s); %Separated by stimulus
                    out = val_stim(:,s);
                    scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),200,out,'s','filled'); axis equal; 
                    title(stim_labels{s});
                end
                
                for r = 1:3
                    subplot(2,3,r+3); %Separated by response
                    out = val_resp(:,r);
                    scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),200,out,'s','filled'); axis equal;
                    title(resp_labels{r});
                end
                
%                 cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
%                     linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
%                 colormap(flipud(cmap));
%                 
                set(get(gcf,'children'),'xtick','','ytick','','box','off','XColor','w','YColor','w');
                set(gcf,'Color','w');
            end
            
            
            %Plot average over mice
            ave_stim = mean(cat(3,val_stim_allSubj{:}),3);
            ave_resp = mean(cat(3,val_resp_allSubj{:}),3);
            figure('name',['Average \Delta' what]);
            
            for s = 1:4
                subplot(2,4,s); %Separated by stimulus
                out = ave_stim(:,s);
                scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),200,out,'s','filled'); axis equal;
                title(stim_labels{s});
            end
            
            for r = 1:3
                subplot(2,3,r+3); %Separated by response
                out = ave_resp(:,r);
                scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),200,out,'s','filled'); axis equal;
                title(resp_labels{r});
            end
            
        end
        
        function plotLogLik(obj)
            %TODO: Change as I now use manual optim for non-laser part
            if isempty(obj.fitData.siteFit)
                error('fit first');
            end
            
            numSubjects = length(obj.names);
            for n = 1:numSubjects
                p_L = obj.fitData.params{n};
                
                figure('name',obj.names{n}); kimg=imread('D:\kirkcaldie_brain_BW.PNG');
                imX=image(-5:1:5,4:-1:-6,kimg); set(gca,'ydir','normal');
                set(gca,'Position',[0 0 1 1]);
                
                [ZL,ZR,~] = obj.getModel(obj.fitData.laserModel,n);
                for site = 1:size(obj.inactivationCoords{n},1)
                    L_idx = obj.data{n}.laserIdx==site;
                    cont = obj.data{n}.contrast_cond(L_idx,:);
                    resp = obj.data{n}.response(L_idx);
                    
                    X = obj.fitData.siteX{site,n};
                    fit = obj.fitData.siteFit{site,n};
                    offset = cvglmnetPredict(fit,X,obj.fitData.lambda,'link'); offset=[offset(:,1)-offset(:,3) offset(:,2)-offset(:,3)];
                    objective = @(PARAMETERS) (obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) );
                    
                    %                     objective([-3000,p_L(site,2),-3000,p_L(site,4)])
                    
                    N = sum(L_idx);
                    fcn_L = @(bL,sL)objective([bL,p_L(site,2),sL,p_L(site,4)])/N ;
                    fcn_R = @(bR,sR)objective([p_L(site,1),bR,p_L(site,3),sR])/N ;
                    
                    % %                     %Plot loglik contour
                    posIn = obj.inactivationCoords{n}(site,1:2)/10 + [0.58 0.465];
                    
                    %                     PLOTRANGES = [0.5*p_L(site,1)-2 0.5*p_L(site,1)+2 0.5*p_L(site,3)-10 0.5*p_L(site,3)+10];
                    PLOTRANGES = [-1 1 -3 3]*4;
                    SAMPLING = 100;
                    ax=axes; hold on;
                    c1=ezcontourf(fcn_L,PLOTRANGES,SAMPLING);
                    
                    if ~isempty(obj.fitData.paramsNull)
                        hold on;
                        scatter(obj.fitData.paramsNull{n}(:,1),obj.fitData.paramsNull{n}(:,3),5,'w.');
                    end
                    
                    set(ax,'Position',[posIn(2)-0.02 posIn(1) 0.05 0.05],'box','off','XColor',[0 0 0 1],'YColor',[0 0 0 1]);
                    set(ax,'color','w');
                    axis square;
                    plot(0,0,'+w');
                    plot(p_L(site,1),p_L(site,3),'ow');
                    hold off;
                    xlabel(''); ylabel('');                     title('');
                    
                    caxis([fcn_L(p_L(site,1),p_L(site,3)) -obj.guess_bpt(n)]);
                    
                    if site == 1
                        xlabel('Bias'); ylabel('Sensitivity'); title('LvNG');
                    else
                        set(gca,'xtick','','ytick','');
                    end
                    
                    %                     PLOTRANGES = [0.5*p_L(site,1)-2 0.5*p_L(site,1)+2 0.5*p_L(site,3)-10 0.5*p_L(site,3)+10];
                    
                    ax=axes; hold on;
                    c2=ezcontourf(fcn_R,PLOTRANGES,SAMPLING);
                    if ~isempty(obj.fitData.paramsNull)
                        hold on;
                        scatter(obj.fitData.paramsNull{n}(:,2),obj.fitData.paramsNull{n}(:,4),5,'w.');
                    end
                    
                    set(ax,'Position',[posIn(2)+0.02 posIn(1) 0.05 0.05],'box','off','XColor',[0 0 0 1],'YColor',[0 0 0 1]);
                    set(ax,'color','w','xtick','','ytick','');
                    axis square;
                    plot(0,0,'+w');
                    plot(p_L(site,2),p_L(site,4),'ow');
                    hold off;
                    xlabel(''); ylabel('');
                    
                    title('');
                    if site == 1
                        title('RvNG');
                    end
                    %                     colormap(flipud(gray));
                    caxis([fcn_R(p_L(site,2),p_L(site,4)) -obj.guess_bpt(n)]);
                    
                    
                    drawnow;
                    
                    %                     keyboard;
                end
            end
        end
        
    end
    
    methods (Access=private)
        function out = cfn(obj,c,subj)
            out = (c.^obj.cfn_parameters{subj}(1))./(c.^obj.cfn_parameters{subj}(1) + obj.cfn_parameters{subj}(2).^obj.cfn_parameters{subj}(1));
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
                load(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' name '.mat']);
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
                                
                                if size(l.inactivationSite,1) ~= 52
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
                
                save(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' name '.mat'],'eRefs');
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
            zl = ZL(offset(:,1),testParams,inputs);
            zr = ZR(offset(:,2),testParams,inputs);
            pL = exp(zl)./(1+exp(zl)+exp(zr));
            pR = exp(zr)./(1+exp(zl)+exp(zr));
            pNG = 1 - pL - pR;
            
            phat = [pL pR pNG];
            
            
        end
        
    end
end