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
    
    properties (Access=public)
        cfn_parameters;
        significance_flag = 0;
        bilateral_flag = 0;
        moreData = 1;
        nestedModel_flag = 0;
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
            %             obj.inactivationCoords = cell(1,numSubjects);
            obj.data = cell(1,numSubjects);
            for n = 1:numSubjects
                D = struct;
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
                        pupil = obj.addPupilData(n,obj.expt);
                        d.pupil = pupil{session};
                    end
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
            
            if ~isempty(strfind(expt,'bilateral'))
                obj.bilateral_flag = 1;
            end
            
            %match laserIdx over all subjects
            %create laserIdx variable
            allD = [obj.data{:}];
            allLas = cat(1,allD.laserCoord);
            [sites,~,ic] = unique(allLas,'rows');
            nanIdx = find(isnan(sites(:,1)),1,'first');
            inactivationSite = sites(1:nanIdx-1,:);
            ic(ic>=nanIdx)=0;
            
            N=arrayfun(@(s)length(s.response),allD);
            for n = 1:numSubjects
                obj.data{n}.laserIdx = ic(1:N(1));
                ic(1:N(1)) = [];
                N(1)=[];
                
                obj.data{n}.areaIdx = obj.identifyArea(obj.data{n}.laserCoord);
            end
            obj.inactivationCoords = inactivationSite;
            
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
            numSubjects = numSubjects + 1;
            
            %now fit c50 and n parameters for each subject's nonlaser data
            obj.cfn_parameters = cell(1,numSubjects);
            for n = 1:numSubjects
                g=GLM(getrow(obj.data{n},obj.data{n}.laserIdx==0)).setModel('C50-subset').fit;
                obj.cfn_parameters{n} = g.parameterFits(5:6);
            end
            obj.plotExptSummary;
        end
        
        function plotExptSummary(obj)
            figure('name','Summary info for data','color','w');
            numSubjects = length(obj.names);
            for n = 1:numSubjects
                tab = tabulate(obj.data{n}.laserIdx);
                subplot(4,numSubjects,n);
                scatter(obj.inactivationCoords(tab(2:end,1),2),obj.inactivationCoords(tab(2:end,1),1),200,tab(2:end,2),'s','filled');
                axis equal;
                title([obj.names{n} ' ' num2str(tab(1,3)) '% noLaser']);
                set(gca,'box','off','xtick','','ytick','','xcolor','w','ycolor','w');
                colorbar;
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
                            numP = 4;
                        case 'sens'
                            ZL = @(off,p,c)(off(:,1) + p(3).*cfn(c(:,1)) + p(5).*cfn(c(:,2)) );
                            ZR = @(off,p,c)(off(:,2) + p(4).*cfn(c(:,2)) + p(6).*cfn(c(:,1)) );
                            numP = 6;
                        case 'sub_sens'
                            ZL = @(off,p,c)(off(:,1) + p(3).*cfn(c(:,1)) );
                            ZR = @(off,p,c)(off(:,2) + p(4).*cfn(c(:,2)) );
                            numP = 4;
                        case 'bias+sub_sens'
                            ZL = @(off,p,c)(off(:,1) + p(1) + p(3).*cfn(c(:,1)) );
                            ZR = @(off,p,c)(off(:,2) + p(2) + p(4).*cfn(c(:,2)) );
                            numP = 4;
                        case 'bias+sens'
                            ZL = @(off,p,c)(off(:,1) + p(1) + p(3).*cfn(c(:,1)) + p(5).*cfn(c(:,2)) );
                            ZR = @(off,p,c)(off(:,2) + p(2) + p(4).*cfn(c(:,2)) + p(6).*cfn(c(:,1)) );
                            numP = 6;

                        case 'nested_bias+sub_sens' %hacky attempt at a nested model implementation
                            ZL = @(off,p,c)( off(:,1) + p(1) + p(3)*max([cfn(c(:,1)) cfn(c(:,2))],[],2)  ); %actually pGO/pNG
                            ZR = @(off,p,c)( off(:,2) + p(2) + p(4)*diff([cfn(c(:,1)) cfn(c(:,2))],[],2) ); %actually pL/pR given GO
                            numP = 4;
                        case 'nested_sub_sens' %hacky attempt at a nested model implementation
                            ZL = @(off,p,c)( off(:,1) + p(3)*max([cfn(c(:,1)) cfn(c(:,2))],[],2)  ); %actually pGO/pNG
                            ZR = @(off,p,c)( off(:,2) + p(4)*diff([cfn(c(:,1)) cfn(c(:,2))],[],2) ); %actually pL/pR given GO
                            numP = 4;
                        case 'nested_bias' %hacky attempt at a nested model implementation
                            ZL = @(off,p,c)( off(:,1) + p(1) ); %actually pGO/pNG
                            ZR = @(off,p,c)( off(:,2) + p(2) ); %actually pL/pR given GO
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
                            numP = 4;
                        case 'sub_sens'
                            ZL = @(p_offset,p,c)( ( p_offset(:,3) + p(3) ).*(cfn(c(:,1)) + p_offset(:,1) ) );
                            ZR = @(p_offset,p,c)( ( p_offset(:,4) + p(4) ).*(cfn(c(:,2)) + p_offset(:,2) ) );
                            numP = 4;
                            
                        case 'nested_bias+sub_sens' %hacky attempt at a nested model implementation
                            ZL = @(p_offset,p,c)( (p_offset(:,3) + p(3)).*(max([cfn(c(:,1)) cfn(c(:,2))],[],2) + p_offset(:,1) + p(1) ) ); %actually pGO/pNG
                            ZR = @(p_offset,p,c)( (p_offset(:,4) + p(4)).*(diff([cfn(c(:,1)) cfn(c(:,2))],[],2) + p_offset(:,2) + p(2)) ); %actually pL/pR given GO
                            numP = 4;
                        case 'nested_bias' %hacky attempt at a nested model implementation
                            ZL = @(p_offset,p,c)( p_offset(:,3).*(max([cfn(c(:,1)) cfn(c(:,2))],[],2) + p_offset(:,1) + p(1) ) ); %actually pGO/pNG
                            ZR = @(p_offset,p,c)( p_offset(:,4).*(diff([cfn(c(:,1)) cfn(c(:,2))],[],2) + p_offset(:,2) + p(2)) ); %actually pL/pR given GO
                            numP = 4;
                        case 'nested_sub_sens' %hacky attempt at a nested model implementation
                            ZL = @(p_offset,p,c)( (p_offset(:,3) + p(3)).*(max([cfn(c(:,1)) cfn(c(:,2))],[],2) + p_offset(:,1) ) ); %actually pGO/pNG
                            ZR = @(p_offset,p,c)( (p_offset(:,4) + p(4)).*(diff([cfn(c(:,1)) cfn(c(:,2))],[],2) + p_offset(:,2)) ); %actually pL/pR given GO
                            numP = 4;
                    end
                    offsetMode = 'P';
                    
            end
            
        end
        
        function obj = fit(obj,laserLambda)
            f=figure;
            obj.fitData.nonLaserModel = 'bias+sub_sens';
            obj.fitData.laserModel = 'bias+sub_sens';
            obj.fitData.biasMode = 'biasAsContrast';
            obj.fitData.laserLambda = laserLambda;
            
            if strcmp(obj.fitData.nonLaserModel(1:4),'nest'); %if nested model, set flag for use later
                obj.nestedModel_flag = 1;
            else
                obj.nestedModel_flag = 0;
            end
            
            numSubjects = length(obj.names);
            obj.fitData.params = cell(1,numSubjects);
            obj.fitData.nonLaserParams = cell(1,numSubjects);
            for n = 1:numSubjects
                
                %First fit non-laser portion of the data
                
                %Manual optim method per-session
                [ZL_nL,ZR_nL,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                nL = obj.data{n}.laserIdx==0;
                
                %                 %Define which parameters are to be regularised
                %                 if numP == 4
                %                                         reg = [3:4];
                % %                     reg = 1:numP;
                %                 else
                %                     reg = 1:numP;
                %                 end
                
                sessions = unique(obj.data{n}.sessionID(nL));
                p_nL = nan(length(sessions),numP);
                for s = 1:length(sessions)
                    cont = obj.data{n}.contrast_cond(nL & obj.data{n}.sessionID==sessions(s),:);
                    resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                    
                    offset = zeros(length(resp),4);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.laserLambda*sum(abs(PARAMETERS)) );
                    
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                            p_nL(s,:) = Opt.solve';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                            [p_nL(s,:),~,exitflag] = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                            
                            if exitflag<1
                                error('Fit error: did not converge');
                                keyboard;
                            end
                    end
                end
                
                obj.fitData.nonLaserParams{n} = p_nL;
                
                %Then go through each inactivation location and fit
                %extended models
                laserIdx = unique(obj.data{n}.laserIdx);
                laserIdx(1) = []; %remove non-laser condition
                [ZL,ZR,numP,offsetMode] = obj.getModel(obj.fitData.laserModel,obj.fitData.biasMode,n);
                p_L = nan(size(obj.inactivationCoords,1),numP); loglik=[];
                for site = unique(laserIdx)'
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
                    objective = @(PARAMETERS) ( -obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) + obj.fitData.laserLambda*sum(abs(PARAMETERS)));
                    
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                            [p,loglik(site,1),~,~] = Opt.solve;
                            p_L(site,:) = p';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                            [p_L(site,:),~,~,~] = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                    end
                    
                    %                     obj.fitData.siteX{site,n} = X;
                    %                     obj.fitData.siteFit{site,n} = fit;
                end
                obj.fitData.params{n} = p_L;
                %                 obj.fitData.paramsLL{n} = -loglik;
            end
            
            disp('Done!'); close(f);
        end
        
        function obj = mergeSites(obj)
            %Merges cortical sites together
            numSubjects = length(obj.names);
                        
            if obj.bilateral_flag==0
                areaList = {'Left visual', [1 2 3 9 10 11 17 18 19];
                            'Right visual', [6 7 8 14 15 16 22 23 24];
                            'Left M2', [42 43 47 48 51];
                            'Right M2', [44 45 49 50 52]};
            else
                areaList = {'Visual',[1 2 3 5 6 7 9 10 11];
                            'M2',[21 22 24 25 26];
                            'Somatosensory',[13 14 15 16 17 18 19 20]};
            end
            
            newCoordList = [];
            for area = 1:size(areaList,1)
                laserIdxs=areaList{area,2};
                coords = obj.inactivationCoords(laserIdxs,:);
%                 warning('INCOMPLETE');
%                 keyboard; 
                for n = 1:numSubjects
                    
                    idx = arrayfun(@(e)any(e==laserIdxs),obj.data{n}.laserIdx);
                    obj.data{n}.laserIdx(idx) = 1000+area;
                end
                
                obj.inactivationCoords(laserIdxs,:)=nan;
                newCoordList = [newCoordList; mean(coords,1)];
            end
            
            %Any remaining laser inactivation trials: merge into 1 misc
            %location
            laserIdxs = find(~isnan(obj.inactivationCoords(:,1)));
            for n = 1:numSubjects
                idx = arrayfun(@(e)any(e==laserIdxs),obj.data{n}.laserIdx);
                obj.data{n}.laserIdx(idx) = 1000+size(areaList,1)+1;
            end
            newCoordList = [newCoordList; 0 0 0];
            obj.inactivationCoords = newCoordList;
            
            %Reset to normal counts
            for n = 1:numSubjects
                obj.data{n}.laserIdx = mod(obj.data{n}.laserIdx,1000);
            end
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
        %
        function obj = fitNull(obj)
            figure;
            kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
            imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
            set(imX,'alphadata',0.7); hold on;
            numSites = size(obj.inactivationCoords,1);
            cols = [linspace(0,1,numSites)' linspace(1,0,numSites)' zeros(numSites,1) ];
            cols(obj.inactivationCoords(:,2)>0,3) = 1;
            scatter(obj.inactivationCoords(:,2),obj.inactivationCoords(:,1),200,cols,'filled');
            ylim([-6 6]);
            
            numSubjects = length(obj.names);
            numIter = 600;
            
            obj.fitData.paramsNull = cell(1,numSubjects);
            figure('name','Sampling from null distribution of laser parameters');
            for n = 1:numSubjects
                %Get all non-laser trials
                tab = tabulate(obj.data{n}.laserIdx);
                proportion = median(tab(2:end,3)/100);
                
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
                    %                     if numP == 4
                    %                         reg = [3:4];
                    %                     else
                    %                         reg = 1:numP;
                    %                     end
                    
                    sessions = unique(obj.data{n}.sessionID(nL));
                    p_nL = nan(length(sessions),numP);
                    for s = 1:length(sessions)
                        cont = obj.data{n}.contrast_cond(nL & obj.data{n}.sessionID==sessions(s),:);
                        resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                        
                        offset = zeros(length(resp),4);
                        objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.laserLambda*sum(PARAMETERS.^2));
                        
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
                                offset(t,:) = [ZL_nL([0 0],p_nL(sess(t),:),cont(t,:)),...
                                    ZR_nL([0 0],p_nL(sess(t),:),cont(t,:)), 0, 0 ];
                            end
                    end
                    
                    [ZL,ZR,numP,offsetMode] = obj.getModel(obj.fitData.laserModel,obj.fitData.biasMode,n);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) + obj.fitData.laserLambda*sum(PARAMETERS.^2));
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
                
                %                 obj.fitData.paramsNullLL{n}=-loglik;
                
                %                 figure; gplotmatrix(obj.fitData.paramsNull{n});
            end
            
            obj.significance_flag = 1;
        end
        
        function obj = removeSite(obj,numSite)
            %Deletes all trials for a particular inactivation site, and
            %adjusts other data structures to ensure no errors occur later
            numSubjects = length(obj.names);
            obj.inactivationCoords(numSite,:) = [];
            for n = 1:numSubjects
                obj.data{n}=getrow(obj.data{n},obj.data{n}.laserIdx~=numSite);
            end
        end
        
        function obj = removeSubj(obj)
            %Deletes all individual subject's data, keeping only pooled data
            numSubjects = length(obj.names);
            obj.expRefs(1:numSubjects-1)=[];
            obj.data(1:numSubjects-1)=[];
            obj.cfn_parameters(1:numSubjects-1)=[];
            obj.guess_bpt(1:numSubjects-1)=[];
            obj.names(1:numSubjects-1)=[];
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
            numSites = size(obj.inactivationCoords,1);
            cols = [linspace(0,1,numSites)' linspace(1,0,numSites)' zeros(numSites,1) ];
            cols(obj.inactivationCoords(:,2)>0,3) = 1;
            %             cols(:,3) = obj.inactivationCoords{1}(:,2)/7 + .5;
            scatter(obj.inactivationCoords(:,2),obj.inactivationCoords(:,1),200,cols,'filled');
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
        
        function crossval_nonLaser(obj)
            %For each mouse, assess different nonlaser models to account
            %for behavioural data (only on nonlaser trials). Here I will
            %test the following models
            %             MNR
            %             Nested logit
            %             C
            %             C50
            %             In each case I will test these models on the entire dataset from each mouse
            numFolds = 10;
            %             models = {'NONLASER_bias+sub_sens_C','NONLASER_bias+sub_sens_C50'};
            models = {'dummy','unconstrained','nested B & C','nested B & C50','mnr B & C','mnr B & C50'};
            numSubjects = length(obj.names);
            figure('name','model comparison','color','w');
            for n = 1:numSubjects
                numSessions = max(obj.data{n}.sessionID);
                
                logLik = nan(numSessions,length(models));
                for s = 1:numSessions
                    D = getrow(obj.data{n},obj.data{n}.laserIdx==0 & obj.data{n}.sessionID==s); %only look at non-laser data
                    D.phat = nan(length(D.response),length(models));
                    [cVal,~,cid] = unique(D.contrast_cond(:));
                    cid = reshape(cid,size(D.contrast_cond));
                    
                    cv = cvpartition(D.response,'KFold',numFolds);
                    
                    for fold = 1:cv.NumTestSets
                        trainC = D.contrast_cond(cv.training(fold),:);
                        trainR = D.response(cv.training(fold));
                        testC = D.contrast_cond(cv.test(fold),:);
                        testR = D.response(cv.test(fold));
                        
                        %                         p_nL = cell(1,length(models));
                        %                         offset = zeros(length(trainR),4);
                        %                         for m = 1:length(models)
                        %                             [ZL,ZR,numP,~] = obj.getModel(models{m},obj.fitData.biasMode,n);
                        %                             objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,trainC,trainR) + obj.fitData.laserLambda*sum(PARAMETERS.^2) );
                        %                             options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000,'Display','off');
                        %                             p_nL{m} = fmincon(objective,zeros(1,numP{1}), [], [], [], [], numP{2}(1,:), numP{2}(2,:), [], options);
                        %
                        %                             %test model
                        %                             p = obj.calculatePhat(ZL,ZR,offset,p_nL{m},testC);
                        % %                             p(p<0.001) = 0.001;
                        %                             D.phat(cv.test(fold),m) = p(:,1).*(testR==1) + p(:,2).*(testR==2) + p(:,3).*(testR==3);
                        %                         end
                        
                        %Test mnr B & C50 model
                        %                         warning('c50 params from full dataset, unfair advantage here');
                        E = struct;
                        E.contrast_cond = trainC;
                        E.response = trainR;
                        g = GLM(E).setModel('C50-subset').fit;
                        cp = g.parameterFits(5:6);
                        cfn = @(c)( (c.^cp(1))./(c.^cp(1) + cp(2).^cp(1)) );
                        
                        %                         cfn = @(c)(obj.cfn(c,n));
                        b = mnrfit(cfn(trainC),trainR);
                        p = mnrval(b,cfn(testC));
                        D.phat(cv.test(fold),6) = p(:,1).*(testR==1) + p(:,2).*(testR==2) + p(:,3).*(testR==3);
                        
                        
                        %Test naive offset only model
                        %                         keyboard;
                        p = sum([trainR==1 trainR==2 trainR==3])/length(trainR);
                        D.phat(cv.test(fold),1) = p(:,1).*(testR==1) + p(:,2).*(testR==2) + p(:,3).*(testR==3);
                        
                        %                         keyboard;
                        %Test nested b & c
                        %                         Z_NG_G = mnrfit(trainC,categorical(trainR<3));
                        Z_NG_G = mnrfit(max(trainC,[],2),categorical(trainR<3));
                        %                         Z_L_R = mnrfit(trainC(trainR<3,:),trainR(trainR<3));
                        Z_L_R = mnrfit(diff(trainC(trainR<3,:),[],2),trainR(trainR<3));
                        %                         pNG = mnrval(Z_NG_G,testC);
                        pNG = mnrval(Z_NG_G,max(testC,[],2));
                        pNG = pNG(:,1);
                        %                         pLR_G = mnrval(Z_L_R,testC);
                        pLR_G = mnrval(Z_L_R,diff(testC,[],2));
                        pLR = bsxfun(@times,1-pNG,pLR_G);
                        p = [pLR pNG];
                        D.phat(cv.test(fold),3) = p(:,1).*(testR==1) + p(:,2).*(testR==2) + p(:,3).*(testR==3);
                        
                        %test nested b & c50
                        xformC = cfn(trainC);
                        %                         Z_NG_G = mnrfit(xformC,categorical(trainR<3));
                        Z_NG_G = mnrfit(max(xformC,[],2),categorical(trainR<3));
                        %                         Z_L_R = mnrfit(xformC(trainR<3,:),trainR(trainR<3));
                        Z_L_R = mnrfit(diff(xformC(trainR<3,:),[],2),trainR(trainR<3));
                        %                         pNG = mnrval(Z_NG_G,cfn(testC)); pNG = pNG(:,1);
                        pNG = mnrval(Z_NG_G,max(cfn(testC),[],2)); pNG = pNG(:,1);
                        
                        %                         pLR_G = mnrval(Z_L_R,cfn(testC));
                        pLR_G = mnrval(Z_L_R,diff(cfn(testC),[],2));
                        pLR = bsxfun(@times,1-pNG,pLR_G);
                        p = [pLR pNG];
                        D.phat(cv.test(fold),4) = p(:,1).*(testR==1) + p(:,2).*(testR==2) + p(:,3).*(testR==3);
                        
                        %Test mnr B & C model
                        b = mnrfit(trainC,trainR);
                        p = mnrval(b,testC);
                        D.phat(cv.test(fold),5) = p(:,1).*(testR==1) + p(:,2).*(testR==2) + p(:,3).*(testR==3);
                        
                        
                        %Test unconstrained model (lookup table of
                        %empirical pL pR pNG values)
                        trainCid = cid(cv.training(fold),:);
                        testCid = cid(cv.test(fold),:);
                        p = nan(length(cVal),length(cVal),3);
                        for cl = 1:length(cVal)
                            for cr = 1:length(cVal)
                                r = trainR(trainCid(:,1)==cl & trainCid(:,2)==cr);
                                p(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                            end
                        end
                        
                        ph = [];
                        for t = 1:length(testR)
                            ph(t,1) = p(testCid(t,1),testCid(t,2),testR(t));
                        end
                        D.phat(cv.test(fold),2) = ph;
                        
%                         if any(isnan(D.phat(cv.test(fold),6)))
%                             keyboard;
%                         end
                    end
                    %                     tab = tabulate(D.response); tab = tab(:,3)/100;
                    %                     guessbpt=sum(tab.*log2(tab));
                    D.phat(D.phat==0)=0.001;
                    logLik(s,:) = nanmean(log2(D.phat),1);
                    
%                     if any(isnan(logLik(s,:)))
%                         keyboard;
%                     end
                end
                subplot(2,round(numSubjects/2),n);
                h=errorbar(1:length(models),mean(logLik),stderr(logLik));
                h.LineStyle='none'; h.Marker='none'; h.LineWidth=3;
                %                 bar(logLik','stacked');
                title(obj.names{n});
                set(gca,'xtick',1:length(models),'XTickLabel',models,'XTickLabelRotation',90,'box','off');
                
                if n == 1
                    ylabel('CV Log_2 likelihood');
                end
                drawnow;
                
%                 keyboard;
            end
        end
        
        function crossval_LaserEffect(obj,n)
            img=imread('D:\kirkcaldie_brain_BW_outline.png');
            numFolds = 10;
            laserModels = {'bias','sub_sens','bias+sub_sens'};
            %             numSubjects = length(obj.names);
            
            figure('name',[obj.names{n} ' ' num2str(numFolds) '-fold CV ' obj.fitData.biasMode],'color','w');
            %             for n = 1:numSubjects
            cv = cvpartition(obj.data{n}.response,'KFold',numFolds);
            phat = nan(length(obj.data{n}.response),length(laserModels)+1);
            for fold = 1:cv.NumTestSets
                
                trainC = obj.data{n}.contrast_cond(cv.training(fold),:);
                trainR = obj.data{n}.response(cv.training(fold));
                trainL = obj.data{n}.laserIdx(cv.training(fold));
                trainSID = obj.data{n}.sessionID(cv.training(fold));
                
                %Fit non-laser model to training set
                [ZL_nL,ZR_nL,numP,offsetMode] = obj.getModel('bias+sub_sens',obj.fitData.biasMode,n);
                nL = (trainL==0);
                
                sessions = unique(trainSID(nL));
                p_nL = nan(length(sessions),numP);
                disp('NON-LASER TRAIN LOGLIK RELATIVE TO GUESS');
                for s = 1:length(sessions)
                    cont = trainC(nL & trainSID==sessions(s),:);
                    resp = trainR(nL & trainSID==sessions(s));
                    
                    offset = zeros(length(resp),4);
                    objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp) + obj.fitData.laserLambda*sum(PARAMETERS.^2) );
                    
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
                p_L = cell(size(obj.inactivationCoords,1),length(laserModels));
                laserIdx = unique(trainL);
                laserIdx(1) = []; %remove non-laser condition
                
                disp('LASER TRAIN LOGLIK RELATIVE TO GUESS');
                for site = laserIdx'
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
            
            phat(phat==0)=0.001;
            
            %                 subplot(3,length(obj.names),n);
            %                 subplot(3,1,1);
            %                 bar(nanmean(log2(phat))-obj.guess_bpt(n));
            %                 set(gca,'XTickLabel',['noLaser' laserModels],'XTickLabelRotation',90,'box','off');
            %                 title(obj.names{n}); ylabel('CV loglik relative to guess');
            
            logLik = nan(size(obj.inactivationCoords,1),size(phat,2));
            laserIdx = unique(obj.data{n}.laserIdx)';
            laserIdx(1) = [];
            for m = 1:size(phat,2)
                for site = laserIdx
                    logLik(site,m) = nanmean(log2(phat(obj.data{n}.laserIdx==site,m)));
                end
                %                     logLik{n}(:,m) = pivottable(obj.data{n}.laserIdx,[],log2(phat(:,m)),'nanmean');
            end
            
            %Plot loglik comparisons
            plotLabels = {'Bias v NoLaser',...
                'Sens v NoLaser',...
                'Bias+Sens v NoLaser',...
                'Bias+Sens v Bias',...
                'Bias+Sens v Sens'};
            plotVal = {logLik(:,2) - logLik(:,1),...
                logLik(:,3) - logLik(:,1),...
                logLik(:,4) - logLik(:,1),...
                logLik(:,4) - logLik(:,2),...
                logLik(:,4) - logLik(:,3)};
            for p = 1:length(plotVal)
                %                     subplot(3,length(obj.names),n+(p*length(obj.names)));
                subplot(round(length(plotVal)/2),2,p); hold on;
                imX=imagesc(linspace(-4.5,4.5,100),linspace(3.75,-5.2,100),img);
                set(gca,'ydir','normal');
                set(imX,'alphadata',0.7);
                
                delta_LL = plotVal{p};
                ap = obj.inactivationCoords(:,1);
                ml = obj.inactivationCoords(:,2);
                if obj.bilateral_flag==1
                    ap = [ap; ap];
                    ml = [ml; -ml];
                    delta_LL = [delta_LL;delta_LL];
                end
                
                h=scatter(ml,ap,200,delta_LL,'s','filled'); caxis([-1 1]*max(abs(cat(1,plotVal{:}))));
                h.MarkerEdgeColor=[1 1 1]*0.95;
                axis equal; set(gca,'xtick','','ytick','','xcolor','w','ycolor','w');
                cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                colormap(flipud(cmap)); colorbar; title(plotLabels{p});
                drawnow;
            end
            
            %             end
        end
        
        function plotFit_NonLaserOneSession(obj,n,sess)
            p_nL = obj.fitData.nonLaserParams{n}(sess,:);
            idx = obj.data{n}.laserIdx==0 & obj.data{n}.sessionID==sess;
            [ZL_nL,ZR_nL,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
            
            cont = obj.data{n}.contrast_cond(idx,:);
            resp = obj.data{n}.response(idx);
            
            cVals = unique(cont(:));
            numPedestals = length(cVals)-1;
            
            figure('name',[obj.names{n} ' session ' num2str(sess)],'color','w');
            %psych curve: 2D representation
            p = nan(length(cVals),length(cVals),3);
            p_pred = nan(size(p));
            for cl = 1:length(cVals)
                for cr = 1:length(cVals)
                    r = resp(cont(:,1) == cVals(cl) & cont(:,2) == cVals(cr));
                    p(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                    p_pred(cl,cr,:) = obj.calculatePhat(ZL_nL,ZR_nL,zeros(1,4),p_nL,[cVals(cl) cVals(cr)]);
                end
            end
            
            labels = {'pL','pR','pNG'};
            for r = 1:3
                subplot(3,4,r);
                imagesc(p(:,:,r)); set(gca,'ydir','normal'); title(['Actual ' labels{r}]);
                caxis([0 1]);
                set(gca,'xtick',1:length(cVals),'ytick',1:length(cVals),'xticklabels',cVals,'yticklabels',cVals);
                xlabel('CR'); ylabel('CL');
                subplot(3,4,r+4);
                imagesc(p_pred(:,:,r)); set(gca,'ydir','normal'); title(['Pred ' labels{r}]);
                caxis([0 1]);
                set(gca,'xtick',1:length(cVals),'ytick',1:length(cVals),'xticklabels',cVals,'yticklabels',cVals);
            end
            axis(get(gcf,'children'),'square');
            %             set(get(gcf,'children'),'ytick',cVals,'xtick',cVals);
            
            
            %pych curve: pedestal representation
            cols = {[0 0.4470 0.7410],...
                [0.8500 0.3250 0.0980],...
                [0.9290 0.6940 0.1250]};
            for ped = 1:numPedestals
                subplot(numPedestals,4,ped*4); hold on;
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
                
                
                ylabel(['ped: ' num2str(cVals(ped))]);
                
                
                if ped ~= numPedestals
                    set(gca,'xtick','','xcolor','w');
                end
                
            end
            
            %chronometry: RTs
            rt = obj.data{n}.RT(idx);
            %             keyboard;
            CL = (cont(:,1) > cont(:,2)) & resp==1;
            CR = (cont(:,1) < cont(:,2)) & resp==2;
            subplot(3,4,9);
            hist(rt(CL)); xlabel('RT left');
            subplot(3,4,10);
            hist(rt(CR)); xlabel('RT right');
            
            set(get(gcf,'children'),'box','off');
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
            
            pCorrect = cell(1,numSubjects);
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
                    
                    % calculate % correct classification in the dataset
                    p_hat = obj.calculatePhat(ZL_nL,ZR_nL,zeros(size(cont,1),4),p_nL,cont);
                    classify = nan(length(resp),1);
                    for t = 1:length(resp)
                        [~,rhat] = max(p_hat(t,:));
                        
                        if rhat == resp(t)
                            classify(t)=1;
                        else
                            classify(t)=0;
                        end
                    end
                    
                    pCorrect{n}(s) = mean(classify);
                    %                     keyboard;
                end
                set(gcf,'color','w');
            end
            
            figure('name','model classification accuracy','color','w');
            for n = 1:numSubjects
                subplot(1,numSubjects,n);
                plot(100*pCorrect{n},'o--'); ylim([0 100]);
                title(obj.names{n});
                set(gca,'box','off');
                
                if n~=1
                    set(gca,'ytick','','ycolor','w');
                else
                    ylabel('% correct'); xlabel('Session');
                end
            end
        end
        
        function plotFit_Laser(obj)
            numSubjects = length(obj.names);
            
            img=imread('D:\kirkcaldie_brain_BW_outline.png');
            
            %Plot laser parameter maps
            figure('name','Laser parameter maps','color','w');
            for n = 1:numSubjects
                p_L = obj.fitData.params{n};
                
                P_old_mag = abs(mean(obj.fitData.nonLaserParams{n},1));
                plotVal = {p_L(:,1),p_L(:,2),p_L(:,3),p_L(:,4)};
%                                 plotVal = {p_L(:,1)/P_old_mag(1),p_L(:,2)/P_old_mag(2),p_L(:,3)/P_old_mag(3),p_L(:,4)/P_old_mag(4)};

                param_labels = {'b_L','b_R','s_L','s_R'};
                
                if obj.significance_flag == 1
                    pvalue = obj.sig{n};
                    dotSize = 35 + (pvalue<0.05)*80 + (pvalue<0.01)*100;
                else
                    dotSize = ones(size(p_L))*200;
                end
                
                for i = 1:length(plotVal)
                    subplot(numSubjects,length(plotVal),length(plotVal)*n-length(plotVal) + i);
                    
                    ap = obj.inactivationCoords(:,1);
                    ml = obj.inactivationCoords(:,2);
                    dot = dotSize(:,i);
                    if obj.bilateral_flag==1
                        ap = [ap; ap];
                        ml = [ml; -ml];
                        plotVal{i} = [plotVal{i}; plotVal{i}];
                        dot = [dot;dot];
                    end
                    
                    imX=imagesc(linspace(-4.5,4.5,100),linspace(3.75,-5.2,100),img);
                    set(gca,'ydir','normal');
                    set(imX,'alphadata',0.7); hold on;
                    
                    
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
                    ylim([-6,4]);
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
            
%             %Plot psych curves showing the laser effect at a few sites, for
%             %each mouse separately
%             
%             if obj.bilateral_flag==0
%                 siteLabels = {'Left V1','Right V1','Left M2','Right M2','Left Barrel','Right Barrel'};
%                 siteID = [10 15 48 49 25 32];
%             elseif obj.bilateral_flag==1
%                 siteLabels = {'V1','M2','Barrel'};
%                 siteID = [7 24 16];
%             end
%             
%             for n = 1:numSubjects
%                 figure('name',[obj.names{n} ' laser psych effects']);
%                 [ZL_nL,ZR_nL,~,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
%                 
%                 for site = 1:length(siteID)
%                     %Plot non-laser curve, use SESSION AVERAGE of parameters for now
%                     p_nL = mean(obj.fitData.nonLaserParams{n});
%                     p_L = obj.fitData.params{n}(siteID(site),:);
%                     
%                     cont = obj.data{n}.contrast_cond;
%                     cVals = unique(cont(:));
%                     numPedestals = length(cVals)-1;
%                     
%                     for ped = 1:numPedestals
%                         subplot(numPedestals,length(siteID),(ped-1)*length(siteID) + site); hold on;
%                         
%                         %Plot non-laser predictions
%                         uC = unique(diff(cont(min(cont,[],2)==cVals(ped),:),[],2));
%                         testCont = [linspace(max(abs(uC))+0.1,0,100)' zeros(100,1); zeros(100,1) linspace(0,max(abs(uC))+0.1,100)'] + cVals(ped);
%                         p_hat = obj.calculatePhat(ZL_nL,ZR_nL,zeros(length(testCont),4),p_nL,testCont);
%                         set(gca,'ColorOrderIndex',1);
%                         plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
%                         
%                         %Try plotting NULL LASER psych bounds
%                         try
%                             p_null = obj.fitData.paramsNull{n};
%                             if size(p_null,2)==2
%                                 p_null = [p_null, zeros(size(p_null))];
%                             end
%                             
%                             p_totalNull = bsxfun(@plus,p_nL,p_null);
%                             numIter = size(p_totalNull,1);
%                             p_hat = nan(size(testCont,1),3,numIter);
%                             for i = 1:numIter
%                                 p_hat(:,:,i) = obj.calculatePhat(ZL_nL,ZR_nL,zeros(length(testCont),4),p_totalNull(i,:),testCont);
%                             end
%                             
%                             ci95 = quantile(p_hat,[0.025 0.975],3);
%                             cols = get(gca,'ColorOrder');
%                             for r = 1:3
%                                 h=fill([diff(testCont,[],2);flipud(diff(testCont,[],2))],[ci95(:,r,1);flipud(ci95(:,r,2))],'b');
%                                 h.EdgeColor=[1 1 1];
%                                 h.FaceAlpha=0.5;
%                                 %                                 h.EdgeAlpha=h.FaceAlpha;
%                                 h.FaceColor=cols(r,:);
%                                 %                                 h.EdgeColor=cols(r,:);
%                             end
%                         catch
%                             warning('Skipping plotting null laser effect bounds');
%                         end
%                         
%                         if site == 1
%                             ylabel(['ped: ' num2str(cVals(ped))]);
%                         else
%                             set(gca,'ytick','','ycolor','w');
%                         end
%                         
%                         if ped ~= numPedestals
%                             set(gca,'xtick','','xcolor','w');
%                         end
%                         
%                         if ped == 1
%                             title(siteLabels{site})
%                         end
%                         
%                         %Overlay laser predictions
%                         if length(p_L)==2
%                             p_L = [p_L 0 0];
%                         end
%                         p_Total = p_nL + p_L;
%                         p_hat = obj.calculatePhat(ZL_nL,ZR_nL,zeros(length(testCont),4),p_Total,testCont);
%                         set(gca,'ColorOrderIndex',1);
%                         plot(diff(testCont,[],2),p_hat,'linewidth',3);
%                         xlim([-1 1]*(max(cVals)+0.1)); ylim([0 1]);
%                     end
%                     
%                 end
%                 
%                 set(gcf,'color','w');
%             end
            
        end
        
        function plotData(obj,what,varargin)
            
            if cell2mat(strfind(varargin,'withSig')) == 1
                withSig = 1;
            else
                withSig = 0;
            end
            
            numSubjects = length(obj.names);
            if cell2mat(strfind(varargin,'pooled')) == 1
                subjList = numSubjects;
            else
                subjList = 1:numSubjects;
            end
            
            img=imread('D:\kirkcaldie_brain_BW_outline.png');
            val_stim_allSubj = cell(1,numSubjects);
            val_resp_allSubj = cell(1,numSubjects);
            
            for n = subjList
                %Get metric of interest for different stimuli and responses
                
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
                stim = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
                stim_labels = {'CL only','CR only','CL = CR','C=[0 0]'};
                
                resp = obj.data{n}.response;
                resp_labels = {'Chose Left','Chose Right','Chose NoGo'};
                laser = obj.data{n}.laserIdx;
                
                switch(what)
                    case 'performance'
                        val = (obj.data{n}.feedbackType==1);
                        cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                            ones(100,1) linspace(1,0,100)' ones(100,1)];
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
                        
                        val = obj.data{n}.RT * 1000;
                        cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                            ones(100,1) linspace(1,0,100)' ones(100,1)];
                    case 'chooseL'
                        val = (resp==1);
                        cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                            ones(100,1) linspace(1,0,100)' ones(100,1)];
                    case 'chooseR'
                        val = (resp==2);
                        cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                            ones(100,1) linspace(1,0,100)' ones(100,1)];
                    case 'chooseNG'
                        val = (resp==3);
                        cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                            ones(100,1) linspace(1,0,100)' ones(100,1)];
                    case 'pupil'
                        val = obj.data{n}.pupil;
                        cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                            linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
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
                val_stim = nan(size(obj.inactivationCoords,1)+1,4);
                val_resp = nan(size(obj.inactivationCoords,1)+1,3);
                for site = 1:(size(obj.inactivationCoords,1)+1)
                    for s = 1:4
                        val_stim(site,s) = nanmean(val(laser==(site-1) & stim==s));
                    end
                    
                    for r = 1:3
                        val_resp(site,r) = nanmean(val(laser==(site-1) & resp==r));
                    end
                end
                
                %Plot map of delta value
                figString = [obj.expt ' ' obj.names{n} ' delta ' what];
                figure('name',figString,'position',[300 500-50*n 1060 550]);
                h=uicontrol('Style','text','String', figString,'Units','normalized','Position', [0 0.9 1 0.1]);
                set(h,'Foregroundcolor','r','FontSize',20,'Fontname','Helvetica','Fontweight','bold');
                
                %compute change in metric from nonLaser condition
                val_stim = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:));
                val_resp = bsxfun(@minus,val_resp(2:end,:),val_resp(1,:));
                
                if withSig == 1
                    %Non-laser sampling to get null distribution of metric
                    NUM_SHUFFLES = 500;
                    val_stim_null = nan(NUM_SHUFFLES,4);
                    val_resp_null = nan(NUM_SHUFFLES,3);
                    tab = tabulate(laser);
                    prop = median(tab(2:end,3)/100); %get avg proportion of laser trials
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
                        
                        for r=1:3
                            a = min([mean(val_resp_null(:,r)>val_resp(t,r)) mean(val_resp_null(:,r)<val_resp(t,r))]);
                            b = min([mean(val_resp_null(:,r)>-val_resp(t,r)) mean(val_resp_null(:,r)<-val_resp(t,r))]);
                            val_resp_pvalue(t,r) = a+b;
                        end
                        
                        for s=1:4
                            a = min([mean(val_stim_null(:,s)>val_stim(t,s)) mean(val_stim_null(:,s)<val_stim(t,s))]);
                            b = min([mean(val_stim_null(:,s)>-val_stim(t,s)) mean(val_stim_null(:,s)<-val_stim(t,s))]);
                            val_stim_pvalue(t,s) = a+b;
                        end
                        
                    end
                    
                else
                    val_resp_pvalue = zeros(size(val_resp,1),3);
                    val_stim_pvalue = zeros(size(val_stim,1),4);
                end
                
                if size(val_stim,2)==3 %if detection task
                    val_stim = [val_stim, zeros(size(val_stim,1),1)];
                end
                
                ap = obj.inactivationCoords(:,1);
                ml = obj.inactivationCoords(:,2);
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
                    imX=imagesc(linspace(-4.5,4.5,100),linspace(3.75,-5.2,100),img); set(gca,'ydir','normal');
                    set(imX,'alphadata',0.7); hold on;
                    
                    h=scatter(ml,ap,dotSize,out,'s','filled'); axis equal;
                    h.MarkerEdgeColor=[1 1 1]*0.95;
                    xlabel(stim_labels{s});
                    %                     caxis([-1 1]*max(abs(val_stim(:,s))));
                    q = max(abs(quantile(val_stim(:),[0.025 0.975])));
                    if isnan(q); q=1; end;
                    caxis([-1 1]*q);
                    colorbar;
                end
                
                for r = 1:3
                    subplot(2,3,r+3); %Separated by response
                    out = val_resp(:,r);
                    dotSize = 25 + (val_resp_pvalue(:,r)<0.05)*140 + (val_resp_pvalue(:,r)<0.01)*165;
                    imX=imagesc(linspace(-4.5,4.5,100),linspace(3.75,-5.2,100),img); set(gca,'ydir','normal');
                    set(imX,'alphadata',0.7); hold on;
                    h=scatter(ml,ap,dotSize,out,'s','filled'); axis equal;
                    h.MarkerEdgeColor=[1 1 1]*0.95;
                    xlabel(resp_labels{r});
                    %                     caxis([-1 1]*max(abs(val_resp(:,r))));
                    q = max(abs(quantile(val_resp(:),[0.025 0.975])));
                    if isnan(q); q=1; end;
                    caxis([-1 1]*q);
                    colorbar;
                end
                
                colormap(cmap);
                %
                set(findobj(get(gcf,'children'),'Type','Axes'),'xtick','','ytick','','box','off','xcolor','k','ycolor','w');
                set(gcf,'Color','w');
                
                val_stim_allSubj{n} = val_stim;
                val_resp_allSubj{n} = val_resp;
                drawnow;
            end
            
            
            %             %Pool data over mice to get average average
            %             val_stim = nan(size(obj.inactivationCoords,1)+1,4);
            %             val_resp = nan(size(obj.inactivationCoords,1)+1,3);
            %             val = cat(1,allSubj.val);
            %             laser = cat(1,allSubj.laser);
            %             resp = cat(1,allSubj.resp);
            %             stim = cat(1,allSubj.stim);
            %
            %             for site = 1:(size(obj.inactivationCoords,1)+1)
            %                 for s = 1:4
            %                     val_stim(site,s) = nanmean(val(laser==(site-1) & stim==s));
            %                 end
            %
            %                 for r = 1:3
            %                     val_resp(site,r) = nanmean(val(laser==(site-1) & resp==r));
            %                 end
            %             end
            %             val_stim = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:));
            %             val_resp = bsxfun(@minus,val_resp(2:end,:),val_resp(1,:));
            %
            %             ap = obj.inactivationCoords(:,1);
            %             ml = obj.inactivationCoords(:,2);
            %             if obj.bilateral_flag==1
            %                 ap = [ap; ap];
            %                 ml = [ml; -ml];
            %                 val_stim = [val_stim; val_stim];
            %                 val_resp = [val_resp; val_resp];
            %             end
            %
            %             %             ave_stim = nanmean(cat(3,val_stim_allSubj{:}),3);
            %             %             ave_resp = nanmean(cat(3,val_resp_allSubj{:}),3);
            %             figString = [obj.expt ' POOLED DATA delta ' what];
            %             figure('name',figString,'position',[1500 200 1060 550]);
            %             h=uicontrol('Style','text','String', figString,'Units','normalized','Position', [0 0.9 1 0.1]);
            %             set(h,'Foregroundcolor','r','FontSize',20,'Fontname','Helvetica','Fontweight','bold');
            %             for s = 1:4
            %                 subplot(2,4,s); %Separated by stimulus
            %                 out = val_stim(:,s);
            %                 scatter(ml,ap,200,out,'s','filled'); axis equal;
            %                 xlabel(stim_labels{s});
            %                 %                 caxis([-1 1]*max(abs(ave_stim(:,s))));
            %                 caxis([-1 1]*max(abs(quantile(val_stim(:),[0.025 0.975]))));
            %                 colorbar;
            %             end
            %
            %             for r = 1:3
            %                 subplot(2,3,r+3); %Separated by response
            %                 out = val_resp(:,r);
            %                 scatter(ml,ap,380,out,'s','filled'); axis equal;
            %                 xlabel(resp_labels{r});
            %                 %                 caxis([-1 1]*max(abs(ave_resp(:,r))));
            %                 caxis([-1 1]*max(abs(quantile(val_resp(:),[0.025 0.975]))));
            %                 colorbar;
            %             end
            %             colormap(cmap);
            %             set(findobj(get(gcf,'children'),'Type','Axes'),'xtick','','ytick','','box','off','xcolor','k','ycolor','w');
            %             set(gcf,'Color','w');
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
                    
                    cVal = unique(D.contrast_cond(:));
                    prop = nan(4,4,3);
                    for cl = 1:length(cVal)
                        for cr = 1:length(cVal)
                            r = D.response(D.contrast_cond(:,1)==cVal(cl) & D.contrast_cond(:,2)==cVal(cr));
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
        
        function crossval_Z(obj,varargin)
            if isempty(varargin)
                DAT = obj.data{end};
                subj = length(obj.names);
            else
                subj = varargin{1};
                DAT = obj.data{subj};
            end
            
            DAT = getrow(DAT,DAT.laserIdx==0);
            DAT.phat = nan(size(DAT.response));
            
            C = cvpartition(DAT.response,'kfold',5);
            for f = 1:C.NumTestSets
                D = getrow(DAT,C.training(f));
                cVal = unique(D.contrast_cond(:));
                p=nan(length(cVal),length(cVal),3);
                for cl = 1:length(cVal)
                    for cr = 1:length(cVal)
                        r = D.response(D.contrast_cond(:,1)==cVal(cl) & D.contrast_cond(:,2)==cVal(cr));
                        p(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                    end
                end
                
                ZL = log(p(:,:,1)./p(:,:,3));
                ZR = log(p(:,:,2)./p(:,:,3));
                
                ZGO = log((1-p(:,:,3))./p(:,:,3));
                ZLRG = log(p(:,:,1)./ p(:,:,2));
                
                idx=reshape(1:16,4,4)';
                ax(1) = subplot(1,2,1);
                ax(2) = subplot(1,2,2); hold(ax(1),'on'); hold(ax(2),'on');
                for i = 1:4
                    plot(ax(1),ZR(idx(i,:)),ZL(idx(i,:)),'ko:');
                    plot(ax(1),ZR(idx(:,i)),ZL(idx(:,i)),'ko:');
                    plot(ax(2),ZLRG(idx(i,:)),ZGO(idx(i,:)),'ko:');
                    plot(ax(2),ZLRG(idx(:,i)),ZGO(idx(:,i)),'ko:');
                end
                
                phat_mnr = nan(size(p));
                phat_mnr(:,:,1) = exp(ZL)./(1+exp(ZL)+exp(ZR));
                phat_mnr(:,:,2) = exp(ZR)./(1+exp(ZL)+exp(ZR));
                phat_mnr(:,:,3) = 1 - phat_mnr(:,:,1) - phat_mnr(:,:,2);
                
                %<= MNR and NESTED REPRESENTATIONS ARE IDENTICAL, BY CONSTRUCTION. TESTING BETWEEN THEM 
                %DOESNT MAKE SENSE ....
%                 pGO = exp(ZGO)./(1+exp(ZGO));
%                 pL = exp(ZLRG)./(1+exp(ZLRG));
%                 pR = 1-pL;
%                 pL = pL.*pGO; pR = pR.*pGO;
%                 pNG = 1-pGO;
%                 phat_nested = nan(size(p));
%                 phat_nested(:,:,1) = pL;
%                 phat_nested(:,:,2) = pR;
%                 phat_nested(:,:,3) = pNG;
                
                %Now test on dataset. Go through each trial of the hold-out
                %set, extract predicted probability of that trial happening
                %based on the training set
                testidx = find(C.test(f));
                D = getrow(DAT,testidx);
                [~,~,cValIdx]=unique(D.contrast_cond);
                cValIdx = reshape(cValIdx,length(D.response),2);
%                 mat = arrayfun(@(r)(phat_mnr(:,:,r)),D.response,'uni',0);
                
                
                for t = 1:length(D.response)
                    DAT.phat(testidx(t)) = phat_mnr(cValIdx(t,1),cValIdx(t,2),D.response(t));
                end
            end
            
            disp(mean(log2(DAT.phat)) - obj.guess_bpt(subj));
        end
        
        function plotData_Z(obj,varargin)
            %Plot empirical Z space with and without laser, for each
            %contrast condition
            
            if isempty(varargin)
                D = obj.data{end};
                subj = length(obj.names);
            else
                subj = varargin{1};
                D = obj.data{subj};
            end
            
            if obj.bilateral_flag==0
                areaList = {
                    'Left visual', [1 2 3 9 10 11 17 18 19]};
%                     'Right visual', [6 7 8 14 15 16 22 23 24];
%                                             'Left M2', [42 43 47 48 51];
%                     'Right M2', [44 45 49 50 52]};
            else
                areaList = {'Visual',[1 2 3 5 6 7 9 10 11];
                    'M2',[21 22 24 25 26]};
%                     'Somatosensory',[13 14 15 16 17 18 19 20]};
            end
            
            if obj.nestedModel_flag == 1
                xlab = 'log(pL/pR) given go';
                ylab = 'log(pGO/pNG)';
            else
                xlab = 'log(pR/pNG)';
                ylab = 'log(pL/pNG)';
            end
            
            
            figure('color','w');
            for area=1:size(areaList,1)
                laserIdxs = areaList{area,2};
                cVal = unique(D.contrast_cond(:));
                
                p=nan(length(cVal),length(cVal),3);
                pLas=nan(length(cVal),length(cVal),3);
                C=nan(length(cVal),length(cVal),2);
                ZL_ci=nan(length(cVal),length(cVal),2);
                ZR_ci=nan(length(cVal),length(cVal),2);
                for cl = 1:length(cVal)
                    for cr = 1:length(cVal)
                        C(cl,cr,:) = cVal([cl cr]);
                        r = D.response(D.contrast_cond(:,1)==cVal(cl) & D.contrast_cond(:,2)==cVal(cr) & D.laserIdx==0);
                        p(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                        
                        psim = squeeze(p(cl,cr,:));
                        n = length(r);
                        psimDraws=mnrnd(n,psim,5000)/n;
                        
                        ZL_ci(cl,cr,:) = quantile(log(psimDraws(:,1)./psimDraws(:,3)),[0.025 0.975]);
                        ZR_ci(cl,cr,:) = quantile(log(psimDraws(:,2)./psimDraws(:,3)),[0.025 0.975]);
                        
                        ZGO_ci(cl,cr,:) = quantile(log((1-psimDraws(:,3))./psimDraws(:,3)),[0.025 0.975]);
                        psimDraws_LR = bsxfun(@rdivide,psimDraws(:,1:2),sum(psimDraws(:,1:2),2));
                        ZLRG_ci(cl,cr,:) = quantile(log(psimDraws_LR(:,1)./psimDraws_LR(:,2)),[0.025 0.975]);
                        
                        r = D.response(D.contrast_cond(:,1)==cVal(cl) & D.contrast_cond(:,2)==cVal(cr) & arrayfun(@(l)(any(l==laserIdxs)),D.laserIdx));
                        pLas(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                        
                    end
                end
                
                
                %                 subplot(3,size(areaList,1),area);
                %                 pL = p(:,:,1);
                %                 pR = p(:,:,2);
                %                 pLLas = pLas(:,:,1);
                %                 pRLas = pLas(:,:,2);
                %                 scatter(pR(:),pL(:),100,'o','filled'); hold on;
                %                 plot(pR([1:4 8 7 6 5 9:12 16 15 14 13]),pL([1:4 8 7 6 5 9:12 16 15 14 13]),'k:');
                %                 plot(pR([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),pL([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),'k:');
                %                 quiver(pR(:),pL(:),pRLas(:)-pR(:),pLLas(:)-pL(:));
                %                 axis equal; xlim([0 1]); ylim([0 1]);
                %                 xlabel('empirical pR'); ylabel('empirical pL'); title(areaList{area,1});
                %
                
                %                 subplot(3,size(areaList,1),area+size(areaList,1)) %Z plot
                h1(area)=subplot(2,size(areaList,1),area);
                
                if obj.nestedModel_flag == 0
                    ZL = log(p(:,:,1)./p(:,:,3));
                    ZR = log(p(:,:,2)./p(:,:,3));
                    
                    ZL_ciLOW = ZL_ci(:,:,1);
                    ZL_ciHIGH = ZL_ci(:,:,2);
                    ZR_ciLOW = ZR_ci(:,:,1);
                    ZR_ciHIGH = ZR_ci(:,:,2);
                    
                    ZLLas = log(pLas(:,:,1)./pLas(:,:,3)); ZLLas(isinf(ZLLas))=NaN;
                    ZRLas = log(pLas(:,:,2)./pLas(:,:,3)); ZRLas(isinf(ZRLas))=NaN;
                    
                    scatter(ZR(:),ZL(:),50,'o','filled'); hold on;
                    plot(ZR([1:4 8 7 6 5 9:12 16 15 14 13]),ZL([1:4 8 7 6 5 9:12 16 15 14 13]),'k:');
                    plot(ZR([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),ZL([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),'k:');
                    quiver(ZR(:),ZL(:),ZRLas(:)-ZR(:),ZLLas(:)-ZL(:),0);
                    %             scatter(ZRLas(:),ZLLas(:),200,'ro','filled'); hold on;
                    %             text(ZR(:),ZL(:),cellfun(@(c)num2str(c),num2cell(1:numSites),'uni',0));
                    ezplot('2*exp(x)=1+exp(x)+exp(y)');
                    ezplot('2*exp(y)=1+exp(x)+exp(y)');
                    ezplot('2=1+exp(x)+exp(y)');
                    
                    l=line([ZR(:)'; ZR(:)'],[ZL_ciLOW(:)';ZL_ciHIGH(:)']); 
                    set(l,'color',[0 0.4470 0.7410],'linestyle','-');
                    l=line([ZR_ciLOW(:)';ZR_ciHIGH(:)'],[ZL(:)'; ZL(:)']); 
                    set(l,'color',[0 0.4470 0.7410],'linestyle','-');
                else
                    pG(:,:,1) = p(:,:,1)./sum(p(:,:,1:2),3);
                    pG(:,:,2) = p(:,:,2)./sum(p(:,:,1:2),3);
                    pGLas(:,:,1) = pLas(:,:,1)./sum(pLas(:,:,1:2),3);
                    pGLas(:,:,2) = pLas(:,:,2)./sum(pLas(:,:,1:2),3);
                    
                    ZGO = log((1-p(:,:,3))./p(:,:,3));
                    ZLRG = log(pG(:,:,1)./pG(:,:,2));
                    
                    ZGO_ciLOW = ZGO_ci(:,:,1);
                    ZGO_ciHIGH = ZGO_ci(:,:,2);
                    ZLRG_ciLOW = ZLRG_ci(:,:,1);
                    ZLRG_ciHIGH = ZLRG_ci(:,:,2);
                    
                    ZGOLas = log((1-pLas(:,:,3))./pLas(:,:,3)); ZGOLas(isinf(ZGOLas))=NaN;
                    ZRLGLas = log(pGLas(:,:,1)./pGLas(:,:,2)); ZRLGLas(isinf(ZRLGLas))=NaN;
                    scatter(ZLRG(:),ZGO(:),50,'o','filled'); hold on;
                    plot(ZLRG([1:4 8 7 6 5 9:12 16 15 14 13]),ZGO([1:4 8 7 6 5 9:12 16 15 14 13]),'k:');
                    plot(ZLRG([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),ZGO([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),'k:');
                    quiver(ZLRG(:),ZGO(:),ZRLGLas(:)-ZLRG(:),ZGOLas(:)-ZGO(:),0);
                    
                    l=line([ZLRG(:)'; ZLRG(:)'],[ZGO_ciLOW(:)';ZGO_ciHIGH(:)']); 
                    set(l,'color',[0 0.4470 0.7410],'linestyle','-');
                    l=line([ZLRG_ciLOW(:)';ZLRG_ciHIGH(:)'],[ZGO(:)'; ZGO(:)']); 
                    set(l,'color',[0 0.4470 0.7410],'linestyle','-');
                    
                    axis equal;
                end
                axis equal;
                xlabel(['empirical ' xlab]); ylabel(['empirical ' ylab]); title(areaList{area,1});
                
                try
                    if cell2mat(strfind(varargin,'withField')) == 1
                        cSet = {cVal,linspace(0,0.54,30)'};
                    else
                        cSet = {cVal};
                    end
                    
                    h2(area)=subplot(2,size(areaList,1),area+size(areaList,1)); hold on;
                    if obj.nestedModel_flag == 0
                        ezplot('2*exp(x)=1+exp(x)+exp(y)');
                        ezplot('2*exp(y)=1+exp(x)+exp(y)');
                        ezplot('2=1+exp(x)+exp(y)');
                    end
                    pNL = mean(obj.fitData.nonLaserParams{subj},1);
                    pL = mean(obj.fitData.params{subj}(laserIdxs,:),1);
                    
                    for i = 1:length(cSet)
                        [cr,cl]=meshgrid(cSet{i});
                        cr = cr(:); cl = cl(:);
                        
                        [ZLfcn,ZRfcn] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,subj);
                        ZL = ZLfcn([0 0 0 0],pNL,[cl cr]); %ZL = reshape(ZL,length(cVal),length(cVal));
                        ZR = ZRfcn([0 0 0 0],pNL,[cl cr]); %ZR = reshape(ZR,length(cVal),length(cVal));
                        
                        [ZLfcn,ZRfcn] = obj.getModel(obj.fitData.laserModel,obj.fitData.biasMode,subj);
                        
                        switch(obj.fitData.biasMode)
                            case 'biasAsOffset' %offset is Z
                                offset = [ZL(:) ZR(:)];
                            case 'biasAsContrast' %offset is params
                                offset = repmat(pNL,length(cl),1);
                        end
                        ZLLas = ZLfcn(offset,pL,[cl cr]); %ZLLas = reshape(ZLLas,length(cVal),length(cVal));
                        ZRLas = ZRfcn(offset,pL,[cl cr]); %ZRLas = reshape(ZRLas,length(cVal),length(cVal));
                        
                        if i ==1 
                            scatter(ZR(:),ZL(:),50,'o','filled');
                            plot(ZR([1:4 8 7 6 5 9:12 16 15 14 13]),ZL([1:4 8 7 6 5 9:12 16 15 14 13]),'k:');
                            plot(ZR([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),ZL([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),'k:');
                            quiver(ZR(:),ZL(:),ZRLas(:)-ZR(:),ZLLas(:)-ZL(:),0);
                        else
                            quiver(ZR(:),ZL(:),ZRLas(:)-ZR(:),ZLLas(:)-ZL(:));
                        end

                    end
                    axis equal;
                    xlabel(['predicted ' xlab]); ylabel(['predicted ' ylab]); title('');
                catch
                end
                
            end
            
            linkaxes([h1 h2],'xy')
            %             end
        end
        
        function plotData_RThist(obj)
            numSubjects = length(obj.names);
            choiceLabel = {'Left choices','Right choices'};
            for n = 1:numSubjects
                for r = 1:2
                    figString = [obj.names{n} ' ' choiceLabel{r}];
                    figure('name',figString,'color','w','WindowScrollWheelFcn',@callback_RTHIST);
                    %                     h=uicontrol('Style','text','String', figString,'Units','normalized','Position', [0 0.9 1 0.1]);
                    %                     set(h,'Foregroundcolor','r','FontSize',20,'Fontname','Helvetica','Fontweight','bold');
                    
                    nonLaserRT = obj.data{n}.RT(obj.data{n}.laserIdx==0 & obj.data{n}.response==r);
                    for site = 1:size(obj.inactivationCoords,1);
                        
                        LaserRT = obj.data{n}.RT(obj.data{n}.laserIdx==site & obj.data{n}.response==r);
                        posIn = obj.inactivationCoords(site,1:2)/10 + [0.43 0.465];
                        axes; hold on;
                        h1 = histogram(log(nonLaserRT),'BinMethod','scott','DisplayStyle','stairs','Normalization','cdf');
                        h2 = histogram(log(LaserRT),'BinMethod','scott','DisplayStyle','stairs','Normalization','cdf');
                        h1.LineWidth=1;
                        h2.LineWidth=1;
                        
                        set(gca,'Position',[posIn(2) posIn(1) 0.07 0.07],'box','off','ytick','');
                        set(gca,'xtick','','ycolor','w','xcolor','w');
                        ylim([0 1]);
                    end
                end
            end
            
            %             allRT = cellfun(@(s)(s.RT),obj.data,'uni',0)'; allRT = cat(1,allRT{:});
            %             allR = cellfun(@(s)(s.response),obj.data,'uni',0)'; allR = cat(1,allR{:});
            %             allL = cellfun(@(s)(s.laserIdx),obj.data,'uni',0)'; allL = cat(1,allL{:});
            %             nonLaserRT = allRT(allL==0 & allR==r);
            %             for r = 1:2
            %                 figString = ['POOLED TRIALS ' choiceLabel{r}];
            %                 figure('name',figString,'color','w','WindowScrollWheelFcn',@callback_RTHIST);
            %
            %                 for site = 1:size(obj.inactivationCoords,1);
            %                     LaserRT = allRT(allL==site & allR==r);
            %                     posIn = obj.inactivationCoords(site,1:2)/10 + [0.43 0.465];
            %                     axes; hold on;
            %                     h1 = histogram(log(nonLaserRT),'BinMethod','scott','DisplayStyle','stairs','Normalization','cdf');
            %                     h2 = histogram(log(LaserRT),'BinMethod','scott','DisplayStyle','stairs','Normalization','cdf');
            %                     h1.LineWidth=1;
            %                     h2.LineWidth=1;
            %
            %                     set(gca,'Position',[posIn(2) posIn(1) 0.07 0.07],'box','off','ytick','');
            %                     set(gca,'xtick','','ycolor','w','xcolor','w');
            %                     ylim([0 1]);
            %                 end
            %             end
            
        end
        
        function plotChoiceEffects(obj)
            %Plot model-free delta pL,pR,pNG, and if available compare with model predicted delta pL,pR,pNG
            cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
                ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
            
            img=imread('D:\kirkcaldie_brain_BW_outline.png');
            
            numSubjects = length(obj.names);
            
            fisherExactP = nan(4,6,2,numSubjects);
            for n = 1:numSubjects
 
                cont = obj.data{n}.contrast_cond;
                c_0 = sum(cont,2)==0;
                c_L = cont(:,1) > 0 & cont(:,2)==0;
                c_R = cont(:,2) > 0 & cont(:,1)==0;
                c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
                stim = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
                stim_labels = {'CL only','CR only','CL=CR','C=[0 0]'};
                col_labels = {'Chose Left','Chose Right','Chose Either'};
                resp = obj.data{n}.response;
                laser = obj.data{n}.laserIdx;
                
%                 %Print Chi2 test for laser effects (collapsing over
%                 %stimuli, area, choices)
%                 areaLabels = {'Left V1','Right V1','Left S1','Right S1','Left M2','Right M2'};
%                 sigString = {'ns','*','**','***','****'};
%                 disp(obj.names{n});
%                 for a = 1:6
%                     disp(['    ' areaLabels{a}]);
%                     for s = 1:4 %test each area
%                         idx = (obj.data{n}.areaIdx==a | obj.data{n}.areaIdx==7); %& (stim==s);
%                         P=[];
%                         [c1,~,P(1)] = crosstab(obj.data{n}.response(idx)==1,obj.data{n}.areaIdx(idx));
%                         [c2,~,P(2)] = crosstab(obj.data{n}.response(idx)==2,obj.data{n}.areaIdx(idx));
% %                         [c3,~,P(3)] = crosstab(obj.data{n}.response(idx),obj.data{n}.areaIdx(idx));
% %                         [c4,~,P(4)] = crosstab(obj.data{n}.response(idx)<3,obj.data{n}.areaIdx(idx));
%                         numStars = 1 + (P(:) < 0.0001) + (P(:) < 0.001) + (P(:) < 0.01) + (P(:) < 0.05);
%                         
%                         [~,P(1),~] = fishertest(c1);
%                         [~,P(2),~] = fishertest(c2);
%                         
%                         disp(['    ' '    ' stim_labels{s}]);
%                         disp(['    ' '    ' '    ' '2x2 (R+NG)/L vs Laser/noLaser : ' sigString{numStars(1)} ' ' num2str(P(1))]);
%                         disp(['    ' '    ' '    ' '2x2 (L+NG)/R vs Laser/noLaser : ' sigString{numStars(2)} ' ' num2str(P(2))]);
%                         %                     disp(['    ' '    ' '3x2   L/R/NG vs Laser/noLaser : ' sigString{numStars(3)} ' ' num2str(P(a,3))]);
%                         %                     disp(['    ' '    ' '2x2 testGO: ' sigString{numStars(4)} ' ' num2str(P(a,4))]);
%                         
%                         fisherExactP(s,a,1,n) = P(1);
%                         fisherExactP(s,a,2,n) = P(2);
%                     end
%                 end
%                 disp(' ');
                
                
                fx1=figure('name',obj.names{n},'color','w');
                resp(stim==0) = [];
                laser(stim==0) = [];
                cont(stim==0,:) = [];
                stim(stim==0)= [];
                
                AX = [];
                %Plot actual data
                val_stim_diffR=[];
                for r = 1:3
                    val = (resp==r);
                    val_stim = nan(size(obj.inactivationCoords,1)+1,4);
                    for site = 1:(size(obj.inactivationCoords,1)+1)
                        for s = 1:4
                            val_stim(site,s) = nanmean(val(laser==(site-1) & stim==s));
                        end
                    end
                    
                    %Change pNoGo to pGo
                    if r==3
                        val_stim = 1 - val_stim;
                    end
                    
                    %compute change in metric from nonLaser condition
                    val_stim_diff = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:));
                    
                    if r==2
                        val_stim_diffR = val_stim_diff;
                    end
                    
                    ap = obj.inactivationCoords(:,1);
                    ml = obj.inactivationCoords(:,2);
                    if obj.bilateral_flag==1
                        ap = [ap; ap];
                        ml = [ml; -ml];
                        val_stim_diff = [val_stim_diff; val_stim_diff];
                    end
                    
                    q = quantile(val_stim_diff(:),[0.025 0.975]);
                    q = max(abs(q));
                    
                    for s = 1:4
                        %                         AX(r,s)=subplot(3,4,s + 4*(r-1));
                        AX(s,r)=subplot(4,3,3*(s-1) + r);
                        hold on;
                        imX=imagesc(linspace(-4.5,4.5,1000),linspace(3.75,-5.2,1000),img,'Parent',AX(s,r));
                        set(gca,'ydir','normal');
                        set(imX,'alphadata',0.7);
                        h=scatter(ml,ap,200,val_stim_diff(:,s),'s','filled'); axis equal;
                        h.MarkerEdgeColor=[1 1 1]*0.95;
                        %                         caxis([-1 1]*q);
                        caxis([-1 1]*0.4);
                        set(gca,'xtick','','ytick','','xcolor','w','ycolor','w');
                        hold on;
                        if r == 1
                            %                             title('Actual & Pred');
                        end
                        
                        if s == 4
                            set(gca,'xcolor','k'); xlabel(col_labels{r});
                        end
                        
                        if r == 1
                            set(gca,'ycolor','k'); ylabel(stim_labels{s});
                        end
                    end
                end
                
                %Plot predictions (ONLY WORKS FOR P OFFSETS CURRENTLY (bias
                %as contrast)
                try
                    nonlaserparams = mean(obj.fitData.nonLaserParams{n},1);
                    laserparams = obj.fitData.params{n};
                    if size(laserparams,2) == 2
                        laserparams = [laserparams zeros(size(laserparams))];
                    end
                    
                    [ZL,ZR,numP,offsetMode] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                    nonLaserP = [];
                    LaserP = [];
                    dR = [];
                    for s = 1:4
                        c = cont(stim==s,:);
                        nonLaserP(:,:,s) = mean(obj.calculatePhat(ZL,ZR,nonlaserparams,[0 0 0 0],c),1);
                        for site = 1:size(obj.inactivationCoords,1)
                            LaserP(site,:,s) = mean(obj.calculatePhat(ZL,ZR,laserparams(site,:)+nonlaserparams,[0 0 0 0],c),1);
                        end
                        
                        %flip pNG to pG
                        nonLaserP(:,3,s) = 1-nonLaserP(:,3,s);
                        LaserP(:,3,s) = 1-LaserP(:,3,s);
                                                
                        dP = bsxfun(@minus,LaserP(:,:,s),nonLaserP(:,:,s));
                        dR = [dR, dP(:,2)];
%                                                 figure(fx2); subplot(1,4,s);
%                                                 plot(val_stim_diff(:,s),dP(:,3),'o'); xlabel('Actual dNG'); ylabel('Pred dNG');
%                                                 set(gca,'box','off');
%                                                 figure(fx1);
%                         
                        ap = obj.inactivationCoords(:,1);
                        ml = obj.inactivationCoords(:,2);
                        if obj.bilateral_flag==1
                            ap = [ap; ap];
                            ml = [ml; -ml];
                            dP = [dP; dP];
                        end
                        
                        for r = 1:3
                            imX=imagesc(linspace(-4.5,4.5,1000)+11,linspace(3.75,-5.2,1000),img,'Parent',AX(s,r)); set(gca,'ydir','normal');
                            set(imX,'alphadata',0.7);
                            h=scatter(AX(s,r),ml+11,ap,200,dP(:,r),'s','filled'); axis equal;
                            h.MarkerEdgeColor=[1 1 1]*0.95;
                        end
                    end
%                     keyboard;
                    %                     titleText = 'Actual & Pred';
                    
                    
                catch
                    
                    %                     titleText = '';
                end
                set(get(gcf,'children'),'box','off');
                colormap(cmap);

                try
                    figure('name',obj.names{n},'color','w');
                    ezplot('y=x'); hold on;
                    
                    plot(val_stim_diffR(:),dR(:),'ko');
                    xlabel('Empirical \Delta p(R)');
                    ylabel('Predicted \Delta p(R)');
                    title('');
                    set(gca,'box','off');
                    xlim([-1 1]); ylim([-1 1]); axis square;
                catch
                end
                
%                Now plot as a scatter plot of pL vs pR
                
                figure('name',obj.names{n},'color','w');
                subplot(4,3,2);
                kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
                imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal','xtick','','ytick','','box','off','xcolor','w','ycolor','w');
                set(imX,'alphadata',0.7); hold on;
                numSites = size(obj.inactivationCoords,1);
                
                %                 cols = get(gca,'ColorOrder');
                %                 cols = [1 0 0;
                %                         0.5 0 0;
                %                         1 0 1;
                %                         0.5 0 0.5;
                %                         0 1 0;
                %                         0 0.5 0;
                %                         0.5 0.5 0.5];
                %                 cols = cols(obj.identifyArea(obj.inactivationCoords),:);
                %                 cols = [linspace(0,1,numSites)' linspace(1,0,numSites)' zeros(numSites,1) ];
                %                 cols(obj.inactivationCoords(:,2)>0,3) = 1;
                %
                
                cols = [69, 198, 234;
                    65, 140, 202;
                    234, 155, 196;
                    176, 115, 175;
                    182, 216, 150;
                    252, 200, 102;
                    201, 54, 0;
                    243, 237, 73;
                    125 125 125]/255;
                
                if obj.bilateral_flag==0
                    areaID = [1 1 2 3 3 2 1 1,...
                        1 1 2 3 3 2 1 1,...
                        5 4 4 3 3 4 4 5,...
                        5 5 5 3 3 5 5 5,...
                        5 5 6 7 7 6 5 5,...
                        6 6 7 7 6 6,...
                        7 7 7 7,...
                        8 8,...
                        9];
                    
                    sym = {'o','o'};
                    
                else
                    areaID = [3 2 1 1,...
                        3 2 1 1,...
                        3 4 4 5,...
                        3 5 5 5,...
                        7 6 5 5,...
                        7 6 6,...
                        7 7,...
                        8];
                    sym = {'o','o'};
                end
                cols = cols(areaID,:);
                
                
                plotCols = cols;
                ap = obj.inactivationCoords(:,1);
                ml = obj.inactivationCoords(:,2);
                Lidx = ml<0;
                Ridx = ~Lidx;
                if obj.bilateral_flag==1
                    ap = [ap; ap];
                    ml = [ml; -ml];
                    plotCols = [plotCols; plotCols];
                end
                
                %                 scatter(ml,ap,200,plotCols,sym{1},'filled');
                h=scatter(ml(Lidx),ap(Lidx),200,plotCols(Lidx,:),sym{1});
                set(h,'linewidth',2);
                scatter(ml(Ridx),ap(Ridx),200,plotCols(Ridx,:),sym{2},'filled');
                ylim([-6 6]);
                %
                val_stim = nan(size(obj.inactivationCoords,1)+1,4,3);
                NL = cell(3,4);
                NLci = cell(3,4);
                for r = 1:3
                    val = (resp==r);
                    for site = 1:(size(obj.inactivationCoords,1)+1)
                        for s = 1:4
                            val_stim(site,s,r) = nanmean(val(laser==(site-1) & stim==s));
                        end
                    end
                    
                    %nonlaser p and pci
                    for s = 1:4
                        nonLasertrials = val(laser==0 & stim==s);
                        [NL{r,s},NLci{r,s}] = binofit(sum(nonLasertrials),length(nonLasertrials));
                    end
                end
                
                for s = 1:4
                    subplot(4,4,s+4); hold on;
                    %                                             scatter(val_stim(2:end,s,2),val_stim(2:end,s,1),50,cols,'filled');
                    h=scatter(val_stim(1+find(Lidx),s,2),val_stim(1+find(Lidx),s,1),50,cols(Lidx,:),sym{1});
                    set(h,'linewidth',2);
                    scatter(val_stim(1+find(Ridx),s,2),val_stim(1+find(Ridx),s,1),50,cols(Ridx,:),sym{2},'filled');
                    
                    scatter(val_stim(1,s,2),val_stim(1,s,1),100,[0 0 0],'filled');
                    %                     l=line(ones(1,2)*NL{2,s},NLci{1,s}); set(l,'linewidth',4,'color',[0 0 0]);
                    %                     l=line(NLci{2,s},ones(1,2)*NL{1,s}); set(l,'linewidth',4,'color',[0 0 0]);
                    
                    pL = val_stim(1,s,1);
                    pR = val_stim(1,s,2);
                    pNG = val_stim(1,s,3);
                    zl = log(pL/pNG);
                    zr = log(pR/pNG);
                    
                    pLn = @(pRn)( 1 - pRn*(1+ exp(-zr)));
                    ezplot(pLn);
                    pLn = @(pRn)( (1 - pRn)/(1+ exp(-zl)));
                    a=ezplot(pLn); a.LineStyle='--';
                    
                    xlabel('pR'); ylabel('pL'); xlim([0 1]); ylim([0 1]); axis square;
                    title(stim_labels{s});
                    
                    try
                        LaserP; %should throw error if model wasn't fitted
                        subplot(4,4,s+8); hold on;
                        h=scatter(LaserP(Lidx,2,s),LaserP(Lidx,1,s),50,cols(Lidx,:),sym{1});
                        set(h,'linewidth',2);
                        scatter(LaserP(Ridx,2,s),LaserP(Ridx,1,s),50,cols(Ridx,:),sym{2},'filled');
                        scatter(nonLaserP(1,2,s),nonLaserP(1,1,s),100,[0 0 0],'filled');
                        pL = nonLaserP(1,1,s);
                        pR = nonLaserP(1,2,s);
                        pNG = nonLaserP(1,3,s);
                        zl = log(pL/pNG);
                        zr = log(pR/pNG);
                        
                        pLn = @(pRn)( 1 - pRn*(1+ exp(-zr)));
                        ezplot(pLn);
                        pLn = @(pRn)( (1 - pRn)/(1+ exp(-zl)));
                        a=ezplot(pLn); a.LineStyle='--';
                        
                        xlabel('pR'); ylabel('pL'); xlim([0 1]); ylim([0 1]); axis square;
                        title('Prediction');
                    catch
                        warning('Need to fit model to plot model predictions');
                    end
                end
                
                %false alarm vs contrtalateral biasing plot
%                 FA_L(n) = val_stim(1,1,1);
%                 FA_R(n) = val_stim(1,1,2);
%                 CL_pR(n) = max(val_stim(:,2,2));
%                 CR_pL(n) = max(val_stim(:,3,1));
                
                %                 %Another view: plot pL vs pR for each possible stimulus
                %                 %configuration, at each site
                % %                 keyboard;
                %                 cVal = unique(cont(:));
                %                 cols = {[3/3 3/3 3/3] [2/3 3/3 2/3] [1/3 3/3 1/3] [0/3 3/3 0/3];
                %                         [3/3 2/3 2/3] [0.78 0.78 0.44] [0.56 0.89 0.22] [1/3 3/3 0/3];
                %                         [3/3 1/3 1/3] [0.89 0.56 0.22] [0.78 0.78 0.11] [2/3 3/3 0/3];
                %                         [3/3 0/3 0/3] [3/3 1/3 0/3] [3/3 2/3 0/3] [3/3 3/3 0/3]};
                %                 [a,b] = meshgrid(cVal);
                %                 val_stim = nan(length(cVal),length(cVal),3,size(obj.inactivationCoords,1)+1);
                %                 for r = 1:3
                %                     val = (resp==r);
                %                     for site = 1:(size(obj.inactivationCoords,1)+1)
                %                         for cl = 1:length(cVal)
                %                             for cr = 1:length(cVal)
                %                                 val_stim(cl,cr,r,site) = nanmean(val(laser==(site-1) & cont(:,1)==cVal(cl) & cont(:,2)==cVal(cr)));
                %                             end
                %                         end
                %                     end
                %                 end
                %
                %                 figure('name',obj.names{n},'color','w');
                % %                 kimg=imread('D:\kirkcaldie_brain_BW.PNG');
                % %                 imX=image(-5:1:5,4:-1:-6,kimg); set(gca,'ydir','normal');
                % %                 set(gca,'Position',[0 0 1 1]);
                %                 for site = 1:size(obj.inactivationCoords,1)
                %                     posIn = obj.inactivationCoords(site,1:2)/10 + [0.58 0.465];
                %                     axes; hold on;
                %
                %                     %plot non-laser dots
                %                     pL = val_stim(:,:,1,1);
                %                     pR = val_stim(:,:,2,1);
                %                     scatter(pR(:),pL(:),50,1-cat(1,cols{:}),'o');
                %
                %                     %Plot site dots
                %                     pL = val_stim(:,:,1,site+1);
                %                     pR = val_stim(:,:,2,site+1);
                %                     scatter(pR(:),pL(:),50,1-cat(1,cols{:}),'o','filled');
                %
                %
                %                     xlabel('pR'); ylabel('pL'); xlim([0 1]); ylim([0 1]); axis square;
                %                     set(gca,'Position',[posIn(2) posIn(1) 0.07 0.07],'box','off','xtick','','ytick','');
                %                 end
                %
                %                 axes; scatter(a(:),b(:),100,1-cat(1,cols{:}),'filled'); xlabel('CR'); ylabel('CL');
                %                 set(gca,'Position',[0.8 0.8 0.07 0.07],'box','off'); axis square; xlim([0 max(cVal)]); ylim([0 max(cVal)]);
            end
            
            %Plot results of significance testing as maps
%             keyboard;
            figure('color','w');
            titles={'(R+NG)/L vs Laser/noLaser','(L+NG)/R vs Laser/noLaser'};
            for n = 1:numSubjects
                
                for i=1:2
                    subplot(numSubjects,2,2*n -2 + i);
                    imagesc((fisherExactP(:,:,i,n)<0.05) + (fisherExactP(:,:,i,n)<0.01) + (fisherExactP(:,:,i,n)<0.001)); axis square;
                    
                    if n<numSubjects
                        set(gca,'xtick','','ytick','');
                    else
                        set(gca,'xtick',1:6,'XTickLabel',areaLabels,'xticklabelrotation',45,'ytick',1:4,'yticklabel',stim_labels);
                    end
                    
                    if n ==1
                        title(titles{i});
                    end
                    
                    if i == 1
                        ylabel(obj.names{n});
                    end
                end
            end
            
            colormap('gray');
            
            %false alarm plots
            %             figure('color','w');
            %             subplot(1,2,1);
            %             plot(FA_L(1:end-1),CR_pL(1:end-1),'o'); xlabel('False alarm left choices'); ylabel('max pL during CR with inactivation');
            %             xlim([0 1]); ylim([0 1]); axis square;
            %             subplot(1,2,2);
            %             plot(FA_R(1:end-1),CL_pR(1:end-1),'o'); xlabel('False alarm right choices'); ylabel('max pR during CL with inactivation');
            %             xlim([0 1]); ylim([0 1]); axis square;
        end
        
        
    end
    
    methods (Access=private)
        function out = cfn(obj,c,subj)
            out = (c.^obj.cfn_parameters{subj}(1))./(c.^obj.cfn_parameters{subj}(1) + obj.cfn_parameters{subj}(2).^obj.cfn_parameters{subj}(1));
        end
        
        function pupil = addPupilData(obj,n,expt)
            
            try
                load(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' obj.names{n} '_PUPIL.mat']);
            catch
                pupil = cell(1,length(obj.expRefs{n}));
                for b = 1:length(obj.expRefs{n})
                    eRef = obj.expRefs{n}{b};
                    block = dat.loadBlock(eRef);
                    trials=block.trial;
                    pupil{b} = nan(block.numCompletedTrials,1);
                    
                    try
                        [mouseName,thisDate,expNum]=dat.parseExpRef(eRef);
                        tsPath = dat.expFilePath(mouseName, thisDate, expNum, 'eyetracking', 'master');
                        load(fullfile(fileparts(tsPath),'eye_processed.mat'));
                        
                        try
                            load(fullfile(fileparts(tsPath), 'eye_timeStamps.mat'));
                        catch
                            
                            %USING NICK'S ALIGNMENT CODE HERE::
                            load(dat.expFilePath(mouseName, thisDate, expNum, 'Timeline', 'master'));
                            tt = Timeline.rawDAQTimestamps;
                            rew = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'rewardEcho'));
                            [~, rewardOnsets] = schmittTimes(tt,rew, [2 3]);
                            aud = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'audioMonitor'));
                            aud = conv(aud.^2, gausswin(100), 'same'); % smooth
                            [~, soundOnsets] = schmittTimes(tt,aud, [0.02 0.03]);
                            % may have to choose this threshold by inspection:
                            % figure; plot(tt, aud);
                            las = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'waveOutput'));
                            [~, laserOnsets] = schmittTimes(tt,aud, [0.2 0.3]);
                            % if using sine wave and only want the first, e.g., then:
                            laserOnsets = laserOnsets(diff([0;laserOnsets])>0.1); % choose events that don't have another event within preceding 100ms
                            alignVideo(mouseName, thisDate, expNum, 'eye');
                            
                            load(fullfile(fileparts(tsPath), 'eye_timeStamps.mat'));
                        end
                        
                        %Go through each trial, and get the eye data at the
                        %right time
                        area = results.area;
                        timestamps = tVid;
                        
                        for t=1:block.numCompletedTrials
                            start = trials(t).onsetToneSoundPlayedTime(1);
                            finish = trials(t).responseMadeTime;
                            idx = (start < timestamps & timestamps < finish);
                            pupil{b}(t) = 2*sqrt(nanmean(area(idx))/pi);
                        end
                        
                        
                        %zscore pupil energy so that pupil diameter is
                        %standardised between sessions
                        pupil{b} = zscore(pupil{b});
                        
                    catch
                        warning('No eye data found.. Setting to NaNs');
                    end
                end
                save(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' obj.names{n} '_PUPIL.mat'],'pupil');
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
                                
                            case 'sparse_unilateral_2D_fullPeriodInactivation'
                                good(b)=0;
                                
                                if length(l.data.response)<150
                                    good(b)=0;
                                end
                                
                                if size(l.inactivationSite,1) < 45
                                    good(b)=0;
                                end
                                
                                if ~any(l.data.laser(:,2)<0)
                                    good(b)=0;
                                end
                                
                                if p.parameters.LaserStimFixedPeriod == 0
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
        
        function phat = calculatePhat(obj,ZL,ZR,offset,testParams,inputs)
            zl = ZL(offset,testParams,inputs);
            zr = ZR(offset,testParams,inputs);
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
        
        function areaID = identifyArea(~,laserCoord)
            %             keyboard;
            areaID = nan(size(laserCoord,1),1);
            for t = 1:length(areaID)
                pos = laserCoord(t,1:2);
                
                if pos(1) <= -2
                    areaID(t) = (pos(2)<0)*1 + (pos(2)>0)*2;
                elseif -1 <= pos(1) && pos(1) <= 1
                    areaID(t) = (pos(2)<0)*3 + (pos(2)>0)*4;
                elseif -2 <= pos(1) && pos(1) <= 3
                    areaID(t) = (pos(2)<0)*5 + (pos(2)>0)*6;
                else
                    areaID(t) = 7;
                end
            end
        end
        
    end
end