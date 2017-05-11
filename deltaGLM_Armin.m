classdef deltaGLM_Armin
    properties
        names;
        data;
        expRefs;
        model;
        perturbation;
        lambda;
        fitData;
        cvData;
        guess_bpt;
        baseDir = '\\basket.cortexlab.net\home\VTA_Inactivation';
        lapseFlag;

    end
    
    properties (Access=public)
    end
    
    methods
        function obj = deltaGLM_Armin
            files = dir([obj.baseDir,'\*.mat']);
            numSubjects = length(files);
            obj.expRefs = cell(1,numSubjects);
            obj.names = cell(1,numSubjects);
            obj.data = cell(1,numSubjects);
            
            for n = 1:numSubjects
                obj.names{n} = files(n).name(1:end-21);
                s = load(fullfile(obj.baseDir,files(n).name));
                obj.expRefs{n} = cellfun(@(str)str(1:end-10),s.list(:),'uni',0);
                
                obj.data{n} = struct;
                for b = 1:length(obj.expRefs{n})
                    block = dat.loadBlock(obj.expRefs{n}{b});
                    trials = block.trial;
                    d = struct;
                    
                    for t=1:block.numCompletedTrials
                        d.stimulus(t,:) = trials(t).condition.visCueContrast';
                        d.response(t,1) = trials(t).responseMadeID';
                        d.repeatNum(t,1) = trials(t).condition.repeatNum;
                        d.reward(t,1) = trials(t).feedbackType==1;
                        d.session(t,1) = b;
                        d.RT(t,1) = trials(t).responseMadeTime - trials(t).interactiveStartedTime;
                        try
                            d.DA(t,1) = trials(t).condition.rewardVolume(2)>0 & trials(t).feedbackType==1;
                        catch
                            d.DA(t,1)= 0;
                        end
                    end
                                        
                    %identify dopamine session type
                    DA_side = unique(d.response(d.DA==1));
                    if isempty(DA_side)
                        DA_side = 0;
                    elseif length(DA_side) == 2
                        DA_side = 3;
                    end
                    d.DAsession = ones(length(d.response),1)*DA_side; %Find trials with DA on the right

                    obj.data{n} = addstruct(obj.data{n},d);
                                        
                    disp([num2str(b) '/' num2str(length(obj.expRefs{n}))]);
                end
                
                %Exclude data where RT > 5sec or repeatnum > 1, or stimulus
                %type is very rare
                [c,ia,ib] = unique(obj.data{n}.stimulus,'rows');
                tab = tabulate(ib);
                stimID_remove = find(tab(:,3)<5);
                idx_remove = any(ib == stimID_remove',2);
                obj.data{n} = getrow(obj.data{n},obj.data{n}.RT < 5 & obj.data{n}.repeatNum == 1 & idx_remove==0);
                
                %compute chance level of classifier
                tab = tabulate(obj.data{n}.response);
                tab = tab(:,3)/100;
                obj.guess_bpt(n)=sum(tab.*log2(tab));
            end

            
            %Create new subject which is the combined data of all subjects
            obj.names = [obj.names '<<ALLSUBJ>>'];
            obj.expRefs = [obj.expRefs {cat(1,obj.expRefs{:})} ];
            obj.data{numSubjects+1} = obj.data{1};
            for n = 2:numSubjects
                addDat = obj.data{n};
                fields = fieldnames(addDat);
                for f = 1:length(fields)
                    obj.data{numSubjects+1}.(fields{f}) = [obj.data{numSubjects+1}.(fields{f}); addDat.(fields{f})];
                end
            end
            tab = tabulate(obj.data{numSubjects + 1}.response);
            tab = tab(:,3)/100;
            obj.guess_bpt(numSubjects + 1)=sum(tab.*log2(tab));
            
        end
        
        function obj = setPerturbation(obj,which)
            
            numSubjects = length(obj.names);
            obj.perturbation = cell(1,numSubjects);
            for n = 1:numSubjects
                
                obj.data{n}.perturbation = nan(size(obj.data{n}.response));
                
                switch(which)
                    case 'DA L/R' %DA on left as 'non perturbation' condition
                        obj.perturbation{n} = {'Left DA sessions','Right DA sessions'};
                        obj.data{n}.perturbation = nan(length(obj.data{n}.response),1);
                        
                        idx = obj.data{n}.DAsession == 1;
                        obj.data{n}.perturbation(idx) = 0;
                        
                        idx = obj.data{n}.DAsession == 2;
                        obj.data{n}.perturbation(idx) = 1;
                        
                    case 'DA R/L' %DA on right as 'non perturbation' condition
                        obj.perturbation{n} = {'Right DA sessions','Left DA sessions'};                        
                        obj.data{n}.perturbation = nan(length(obj.data{n}.response),1);
                        
                        idx = obj.data{n}.DAsession == 1;
                        obj.data{n}.perturbation(idx) = 1;
                        
                        idx = obj.data{n}.DAsession == 2;
                        obj.data{n}.perturbation(idx) = 0;
                    case 'DA L/R per session'
                        %No perturbation: DA_L sessions all grouped
                        %together
                        %Perturbation: DA_R sessions fit individually
                        obj.data{n}.perturbation = nan(length(obj.data{n}.response),1);
                        
                        idx = obj.data{n}.DAsession == 1;
                        obj.data{n}.perturbation(idx) = 0;
                        
                        idx = obj.data{n}.DAsession == 2;
                        sess = unique(obj.data{n}.session(idx));
                        
                        for s = 1:length(sess)
                            idx = obj.data{n}.DAsession == 2 & obj.data{n}.session == sess(s);
                            obj.data{n}.perturbation(idx) = s;
                        end
                        
                        labels = arrayfun(@(id)['Right DA ' num2str(id)],sess,'uni',0);
                        obj.perturbation{n} = ['Left DA sessions',labels'];
                        
                    otherwise
                        error('Doesnt exist');
                end

            end
            
        end

        
        function m = getModel(~,model,deltaStr)
           m.name = {model,deltaStr};
            switch(model)
                case '2AFC-BO'
                    m.Z = @(p_offset,p,c)( p(1) );
                    m.LB = [-inf ]; %only used when doing non-laser fit
                    m.UB = [+inf ]; %only used when doing non-laser fit
                    m.pLabels = {'b'};
                    
                    m.deltaZ = @(p_offset,p,c)( p_offset(:,1) + p(1) );
                    
                case '2AFC-BO_C'
                    m.Z = @(p_offset,p,c)( p(1) + 10*p(2).*c(:,1) + 10*p(3).*c(:,2) );
                    m.LB = [-inf -inf -inf]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf]; %only used when doing non-laser fit
                    m.pLabels = {'b','sL','sR'};
                    
                    m.deltaZ = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,2)*2^p(2)).*c(:,1) + 10*(p_offset(:,3)*2^p(3)).*c(:,2)   );
                    
                case '2AFC-BO_C^N'
                    m.cfn = @(c,n)(c.^n);
                    m.Z = @(p_offset,p,c)( p(1) + 10*p(2).*m.cfn(c(:,1),p(4)) + 10*p(3).*m.cfn(c(:,2),p(4)) );
                    m.LB = [-inf -inf -inf 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf 1 0]; %only used when doing non-laser fit
                    m.pLabels = {'b','sL','sR','nL','nR'};
                    
                    m.deltaZ = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,2)*2^p(2)).*m.cfn(c(:,1),p_offset(:,4)*2^p(4)) + 10*(p_offset(:,3)*2^p(3)).*m.cfn(c(:,2),p_offset(:,5)*2^p(5))     );
                    
                case '2AFC-BO_NC50'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.Z = @(p_offset,p,c)( p(1) + 10*p(2).*m.cfn(c(:,1),p(4),p(6)) + 10*p(3).*m.cfn(c(:,2),p(4),p(6)) );
                    m.LB = [-inf -inf    -inf    0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf   20  0 3     0]; %only used when doing non-laser fit
                    m.pLabels = {'b','sL','sR','nL','nR','c50L','c50R'};
                    
                    m.deltaZ = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,2)*2^p(2)).*m.cfn(c(:,1),p_offset(:,4)*2^p(4),p_offset(:,6)*2^p(6)) + 10*(p_offset(:,3)*2^p(3)).*m.cfn(c(:,2),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) );
                    
                case '2AFC-BO_NC50DIFF'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.Z = @(p_offset,p,c)( p(1) + 10*p(2).*(m.cfn(c(:,1),p(3),p(5)) - m.cfn(c(:,2),p(3),p(5))) );
                    m.LB = [-inf -inf 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf 20  0 3     0]; %only used when doing non-laser fit
                    m.pLabels = {'b','s','nL','nR','c50L','c50R'};
                    
                    m.deltaZ = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,2)*2^p(2)).*(m.cfn(c(:,1),p_offset(:,3)*2^p(3),p_offset(:,5)*2^p(5)) - m.cfn(c(:,2),p_offset(:,4)*2^p(4),p_offset(:,6)*2^p(6))) );
                    
                otherwise
                    error(['Model ' model ' does not exist']);
            end
            
            if isempty(deltaStr)
                m.deltaP = zeros(1,length(m.pLabels));
            else
                m.deltaP = zeros(1,length(m.pLabels));

                if ischar(deltaStr)
                    if any(strfind(deltaStr,'bias'))
                        m.deltaP(1) = 1;
                    end
                    
                    if any(strfind(deltaStr,'sens'))
                        m.deltaP(2:3) = 1;
                    end
                    
                    if any(strfind(deltaStr,'shape')) 
                        m.deltaP(4:end) = 1;
                    end
                    
                    if any(strfind(deltaStr,'^n'))
                        m.deltaP(4:5) = 1;
                    end
                    
                    if any(strfind(deltaStr,'^c50'))
                        m.deltaP(6:7) = 1;
                    end
                    
                    if any(strfind(deltaStr,'all'))
                        m.deltaP(1:end) = 1;
                    end
                else
                    m.deltaP = deltaStr;
                end
            end
        end

        function obj = fit(obj,behavModel,perturbationModel,withLapse)
            obj.model = obj.getModel(behavModel,perturbationModel);

            if withLapse == 1
                obj.model.pLabels = [obj.model.pLabels 'lapse'];
                obj.model.LB = [obj.model.LB, 0];
                obj.model.UB = [obj.model.UB, 1];
                
                obj.lapseFlag = 1;
                
                if any(strfind(perturbationModel,'lapse')) 
                   obj.model.deltaP = [obj.model.deltaP 1];
                else
                   obj.model.deltaP = [obj.model.deltaP 0]; 
                end
            else
                obj.lapseFlag = 0;
            end

            obj.lambda = 0.001;
                       
            numSubjects = length(obj.names);
            obj.fitData.params = cell(1,numSubjects);
            obj.fitData.deltaParams = cell(1,numSubjects);
            for n = 1:numSubjects

                idx = obj.data{n}.perturbation == 0;
                cont = obj.data{n}.stimulus(idx,:);
                resp = obj.data{n}.response(idx);
                p_nL = obj.util_optim(obj.model.Z,[],cont,resp,obj.model.LB,obj.model.UB);
                
                p_nL(:,obj.model.LB==0 & obj.model.UB==0) = p_nL(:,find(obj.model.LB==0 & obj.model.UB==0) - 1);
                
                obj.fitData.params{n} = p_nL;
                
                %Then go through each perturbation condition and fit
                %extended models
                
                %Force UB and LB to be 0 for parameters which are NOT
                %changed by the perturbations
                LB = -inf(1,length(obj.model.pLabels));
                UB = inf(1,length(obj.model.pLabels));
                LB(~obj.model.deltaP)=0;
                UB(~obj.model.deltaP)=0;
                
                perturbTrials = obj.data{n}.perturbation~=0 & ~isnan(obj.data{n}.perturbation);
                perturbConds = unique(obj.data{n}.perturbation(perturbTrials));
                
                p_L = [];
                for p = 1:length(perturbConds)
                    idx = obj.data{n}.perturbation==perturbConds(p);
                    cont = obj.data{n}.stimulus(idx,:);
                    resp = obj.data{n}.response(idx);
                    offset = p_nL;
                    p_L(p,:) = obj.util_optim(obj.model.deltaZ,offset,cont,resp,LB,UB);
                end

                obj.fitData.deltaParams{n} = p_L;                
            end

            disp('Done!'); 
            
            obj.plotFit_mega;
        end
        
        function plotFit_mega(obj)
            numSubjects = length(obj.names)-1;
            numParams = length(obj.model.pLabels);
            
            figure('name','delta parameters','color','w');
            
            for p = 1:numParams
                subplot(1,numParams,p);
                title(obj.model.pLabels{p});
                hold on;
                for n = 1:numSubjects
                    dP = obj.fitData.deltaParams{n}(:,p);
                    h = notBoxPlot(dP,ones(length(dP),1)*n);
                    h.data.Marker='.';
                    h.data.MarkerSize=10;
                end
                set(gca,'xtick',1:numSubjects,'xticklabels',obj.names(1:end-1),'xticklabelrotation',45);
                
                l= line([0 numSubjects],[0 0]);
                l.LineWidth=1;
                hold off;
                
            end
            
        end
        
        function likelihoodratiotest(obj,baseModel,pmodel_unrestricted,pmodel_restricted)
            obj.lambda = 0;
            
            base = obj.getModel(baseModel,[]);
            pModel{1} = obj.getModel(baseModel,pmodel_unrestricted);
            pModel{2} = obj.getModel(baseModel,pmodel_restricted);
            
            
            numSubjects = length(obj.names);
            for n = 1:numSubjects
                disp(obj.names{n});
                %First fit base non-perturbation portion of the data
                
                %Manual optim method per-session
                nL = obj.data{n}.perturbation==0;
                
                cont = obj.data{n}.stimulus(nL,:);
                resp = obj.data{n}.response(nL);
                
                p_nL = obj.util_optim(base.Z,[],cont,resp,base.LB,base.UB);
                
                
                %If the base model does not fit shape params on each
                %side, make sure to copy the parameter used for each side
                p_nL(base.LB==0 & base.UB==0) = p_nL(find(base.LB==0 & base.UB==0) - 1);
                
                
                %Then go through each perturbation condition and fit
                %extended models
                
                pCond = (1:length(obj.perturbation))-1; pCond(1)=[];
                
                for o = 1:length(pCond)
                    L_idx = obj.data{n}.perturbation==pCond(o);
                    
                    cont = obj.data{n}.stimulus(L_idx,:);
                    resp = obj.data{n}.response(L_idx);
                    offset = repmat(p_nL,length(resp),1);
                    
                    loglik=[];
                    for id = 1:2
                        LB = -inf(1,length(pModel{id}.pLabels));
                        UB = inf(1,length(pModel{id}.pLabels));
                        LB(~pModel{id}.deltaP)=0;
                        UB(~pModel{id}.deltaP)=0;
                        
                        p_L = obj.util_optim(pModel{id}.deltaZ,offset,cont,resp,LB,UB);
                        phat = obj.calculatePhat(pModel{id}.deltaZ,offset,p_L,cont);
                        p = (resp==1).*phat(:,1) + (resp==2).*phat(:,2);
                        
                        loglik(id) = sum(log(p));
                    end
                    
                    dof = sum(pModel{1}.deltaP) - sum(pModel{2}.deltaP);
                    if dof < 0
                        error('Badly ordered arguments. First: unrestricted model, second: restricted model');
                    end
                    
                    [h,pValue,stat,cValue] = lratiotest(loglik(1),loglik(2),dof,0.05);
                    if pValue > 0.05
                        symbol = '';
                    elseif pValue < 0.001
                        symbol = '***';
                    elseif pValue < 0.01
                        symbol = '**';
                    elseif pValue < 0.05
                        symbol = '*';
                    end
                    disp([obj.perturbation{pCond(o)} '. LogLiks:[' num2str(loglik(1)) ',' num2str(loglik(2)) '] ' symbol]);
                end
                
                disp(' ');
            end
      
            
            
        end

        function plotFit(obj)
            numSubjects = length(obj.names);
            
            %PLOT NON-PERTURBATION FIT + DATA
            cols = {[0 0.4470 0.7410],...
                [0.8500 0.3250 0.0980],...
                [0.9290 0.6940 0.1250]};
            
            for n = 1:numSubjects
                figure('name',[obj.names{n} ' predictions against actual psychometric data'],'color','w');
                idx = obj.data{n}.perturbation == 0;
                
                cont = obj.data{n}.stimulus(idx,:);
                resp = obj.data{n}.response(idx);
                
                cVals = unique(cont,'rows');
                
                perturbTrials = obj.data{n}.perturbation~=0 & ~isnan(obj.data{n}.perturbation);
                perturbConds = unique(obj.data{n}.perturbation(perturbTrials));
                
                numDeltaConditions = size(perturbConds,1);
                subplot(2,numDeltaConditions+1,1);
                hold on; title(obj.perturbation{n}{1});
                
                ph=[];
                for c = 1:size(cVals,1)
                    r = resp(cont(:,1)==cVals(c,1) & cont(:,2)==cVals(c,2));
                    [ph(c,:),pci] = binofit(sum([r==1 r==2],1),length(r));
                    
                    dC = cVals(c,2)-cVals(c,1);
                    for ch=1:2
                        l(1)=line([1 1]*dC,pci(ch,:));
                        %                             l(2)=line([dC-0.03 dC+0.03],[1 1]*pci(ch,1));
                        %                             l(3)=line([dC-0.03 dC+0.03],[1 1]*pci(ch,2));
                        set(l,'Color',cols{ch},'Linewidth',0.5);
                    end
                end
                plot(cVals(:,2)-cVals(:,1),ph,'.','markersize',15);
                
                %Plot predictions
                testCont = [linspace(1,0,100)' zeros(100,1); zeros(100,1) linspace(0,1,100)'];
                p_hat = obj.calculatePhat(obj.model.Z,[],obj.fitData.params{n},testCont);
                set(gca,'ColorOrderIndex',1);
                plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
                ylim([0 1]);
                
                
                
                %PLOT PERTURBATION DATASETS
                for p = 1:length(perturbConds)
                    subplot(2,numDeltaConditions+1,1+p);
                    hold on;
                    
                    idx = obj.data{n}.perturbation == perturbConds(p);
                    cont = obj.data{n}.stimulus(idx,:);
                    resp = obj.data{n}.response(idx);
                    cVals = unique(cont,'rows');
                    
                    ph=[];
                    for c = 1:size(cVals,1)
                        r = resp(cont(:,1)==cVals(c,1) & cont(:,2)==cVals(c,2));
                        [ph(c,:),pci] = binofit(sum([r==1 r==2],1),length(r));
                        
                        dC = cVals(c,2)-cVals(c,1);
                        for ch=1:2
                            l(1)=line([1 1]*dC,pci(ch,:));
                            set(l,'Color',cols{ch},'Linewidth',0.5);
                        end
                    end
                    plot(cVals(:,2)-cVals(:,1),ph,'.','markersize',15);
                    
                    %Plot predictions
                    testCont = [linspace(1,0,100)' zeros(100,1); zeros(100,1) linspace(0,1,100)'];
                    p_hat = obj.calculatePhat(obj.model.deltaZ,obj.fitData.params{n},obj.fitData.deltaParams{n}(p,:),testCont);
                    set(gca,'ColorOrderIndex',1);
                    plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
                    ylim([0 1]);
                    
                    title(obj.perturbation{n}{p+1});
                    
                    %PLOT PARAMETER PERTURBATIONS
                    subplot(2,numDeltaConditions+1,1+p+numDeltaConditions+1);
                    h=bar(obj.fitData.deltaParams{n}(p,:));
                    h.EdgeColor='none';
                    set(gca,'XTickLabel',obj.model.pLabels,'XTickLabelRotation',45);
                    
                    
                end

            end
            
            
            
            
        end
 
        function crossval_basemodel(obj)
            
            numSubjects = length(obj.names);
            for n = 1:numSubjects
                
                try
                    
                    D = getrow(obj.data{n},obj.data{n}.DAsession==0);
                    if isempty(D.response)
                        error('No non-DA sessions exist for this subject');
                    end
                    
                    sessions = unique(D.session);
                    
                    %For each mouse, assess different nonlaser models to account
                    %for behavioural data (only on nonlaser trials).
                    
                    
                    numFolds = 5;
                    models = {'2AFC-BO','2AFC-BO_C','2AFC-BO_C^N','2AFC-BO_NC50'};
                    
                    %             cfn = @(c,n,c50) ( (c.^n)./(c50.^n + c.^n) );
                    
                    %             axes; hold on;
                    
                    logLik_perSession = nan(length(sessions),length(models));
                    for session=1:length(sessions)
                        disp([num2str(session) '/' num2str(length(sessions))]);
                        E = getrow(D,D.session==sessions(session));
                        E.phat = nan(length(D.response),length(models));
                        
                        cv=cvpartition(E.response,'kfold',numFolds);
                        for fold = 1:cv.NumTestSets
                            trainC = E.stimulus(cv.training(fold),1:2);
                            trainR = E.response(cv.training(fold));
                            testC = E.stimulus(cv.test(fold),1:2);
                            testR = E.response(cv.test(fold));
                            
                            for m = 1:length(models)
                                %                         disp(m);
                                
                                model = obj.getModel(models{m},[]);
                                
                                [PARAM,exitflag] = obj.util_optim(model.Z,[],trainC,trainR,model.LB,model.UB);
                                PARAM(:,model.LB==0 & model.UB==0) = PARAM(:,find(model.LB==0 & model.UB==0) - 1);
                                if exitflag == 1 %only test if model succeeded in fitting
                                    p = obj.calculatePhat(model.Z,[],PARAM,testC);
                                    E.phat(cv.test(fold),m) = p(:,1).*(testR==1) + p(:,2).*(testR==2);
                                end
                            end
                        end
                        
                        logLik_perSession(session,:) = nanmean(log2(E.phat),1);
                    end
                    figure('name','model comparison','color','w'); 
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
                    
                    
                    
                    
                catch
                end
            end
            
            
        end
        
    end
    
    methods (Access=private)
     
        function [p,exitflag] = util_optim(obj,Z,offset,cont,resp,LB,UB)

            objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,Z,offset,cont,resp) + obj.lambda*sum(abs(PARAMETERS).^2) );
            
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
%                 keyboard;
            end
        end
        
        function logLik = calculateLogLik(obj,testParams,Z,offset,inputs,responses)
            phat = obj.calculatePhat(Z,offset,testParams,inputs);
            p = (responses==1).*phat(:,1) + (responses==2).*phat(:,2);
            
            p(p<0.0001) = 0.0001;
            logLik = nanmean(log2(p));

        end
        
        function phat = calculatePhat(obj,Z,offset,testParams,inputs)
            switch(obj.lapseFlag)
                case 0
                    lapse = 0;
                case 1
                    lapse = testParams(end);
            end
            
            z = Z(offset,testParams,inputs);
            pL = lapse + (1-2*lapse)*exp(z)./(1+exp(z));
            pR = 1-pL;
            
            phat = [pL pR];
        end
        
        
    end
end