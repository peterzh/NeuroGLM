classdef deltaGLM
    properties
        names;
        data;
        expRefs;
        model;
        lambda;
        fitData;
        cvData;
        guess_bpt;
        baseDir = '\\zserver.cortexlab.net\Lab\Share\Miles\adaptation';
    end
    
    properties (Access=public)
    end
    
    methods
        function obj = deltaGLM(names)
            numSubjects = length(names);
            obj.data = cell(1,numSubjects);
            obj.names = names;
            obj.expRefs = cell(1,numSubjects);
            
            for n = 1:numSubjects
                
                d = load(fullfile(obj.baseDir,[names{n} '-data-all.mat']));
                d = d.(names{n});
                
                D = struct;
                D.stimulus = d.stimulus;
                D.stimulusOrientation = d.testOrient;
                D.response = d.response;
                D.repeatNum = d.repeatNum;
%                 D.RT = d.rt';
                D.adapterContrast = d.adapterContrast;
                D.adapterOrientation = d.adapterOrient;
                D.adapterCongruence = D.adapterOrientation - D.stimulusOrientation;
                
                %If response is in [-1 0 1], convert to [1 3 2]
                if min(D.response)==-1
                    D.response(D.response==0)=3;
                    D.response(D.response==1)=2;
                    D.response(D.response==-1)=1;
                end
                
                
                obj.data{n} = D;
                obj.expRefs{n} = d.expRef;
                
                tab = tabulate(D.response);
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

        
        function m = getModel(~,model,deltaStr)
           m.name = {model,deltaStr};
            switch(model)
                case 'BO_C^N'
                    m.cfn = @(c,n)(c.^n);
                    m.ZL = @(p_offset,p,c)( p(1) + 10*p(3).*m.cfn(c(:,1),p(5)) );
                    m.ZR = @(p_offset,p,c)( p(2) + 10*p(4).*m.cfn(c(:,2),p(5)) );
                    m.LB = [-inf -inf 0 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 1 0]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','sL','sR','nL','nR'};
                    
                    m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,3)*2^p(3)).*m.cfn(c(:,1),p_offset(:,5)*2^p(5)) );
                    m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + 10*(p_offset(:,4)*2^p(4)).*m.cfn(c(:,2),p_offset(:,6)*2^p(6))  );
                    
                case 'BC_C^N'
                    m.cfn = @(c,n)(c.^n);
                    m.ZL = @(p_offset,p,c)( 10*p(3).*(m.cfn(c(:,1),p(5)) + p(1)) );
                    m.ZR = @(p_offset,p,c)( 10*p(4).*(m.cfn(c(:,2),p(5)) + p(2)) );
                    m.LB = [-inf -inf 0 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 1 0]; %only used when doing non-laser fit
                    m.pLabels = {'bL','bR','sL','sR','nL','nR'};
                    
                    m.deltaZL = @(p_offset,p,c)( 10*(p_offset(:,3)*2^p(3)).*(m.cfn(c(:,1),p_offset(:,5)*2^p(5)) + p_offset(:,1) + p(1)) );
                    m.deltaZR = @(p_offset,p,c)( 10*(p_offset(:,4)*2^p(4)).*(m.cfn(c(:,2),p_offset(:,6)*2^p(6)) + p_offset(:,2) + p(2)) );
                    
                case 'BO_NC50'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p_offset,p,c)( p(1) + 10*p(3).*m.cfn(c(:,1),p(5),p(7)) );
                    m.ZR = @(p_offset,p,c)( p(2) + 10*p(4).*m.cfn(c(:,2),p(5),p(7)) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'bL','bR','sL','sR','nL','nR','c50L','c50R'};

                    m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,3)*2^p(3)).*m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) );
                    m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + 10*(p_offset(:,4)*2^p(4)).*m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8))  );

                case 'BC_NC50'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
%                     m.cfn = @(c,varargin)((c.^varargin{1})./(c.^varargin{1} + varargin{2}.^varargin{1}));
                    m.ZL = @(p_offset,p,c)( 10*p(3).*(m.cfn(c(:,1),p(5),p(7)) + p(1)) );
                    m.ZR = @(p_offset,p,c)( 10*p(4).*(m.cfn(c(:,2),p(5),p(7)) + p(2)) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'bL','bR','sL','sR','nL','nR','c50L','c50R'};
                    
                    m.deltaZL = @(p_offset,p,c)( 10*(p_offset(:,3)*2^p(3)).*(m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) + p_offset(:,1) + p(1)) );
                    m.deltaZR = @(p_offset,p,c)( 10*(p_offset(:,4)*2^p(4)).*(m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)) + p_offset(:,2) + p(2)) );
                    
                case 'BO_NC50_NESTED'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p_offset,p,c)( p(1) + p(3).*max([m.cfn(c(:,1),p(5),p(7)) m.cfn(c(:,2),p(5),p(7))],[],2) ); %actually ZGnG
                    m.ZR = @(p_offset,p,c)( p(2) + p(4).*diff([m.cfn(c(:,1),p(5),p(7)) m.cfn(c(:,2),p(5),p(7))],[],2) ); %actually ZLR
                    m.LB = [-inf -inf -inf -inf 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'bGO','bLR','sGO','sLR','nL','nR','c50L','c50R'};
                    
                    m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + (p_offset(:,3)*2^p(3)).*max([m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) m.cfn(c(:,2),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7))],[],2) ); %actually ZGnG
                    m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + (p_offset(:,4)*2^p(4)).*diff([m.cfn(c(:,1),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)) m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8))],[],2) ); %actually ZLR

                otherwise
                    error(['Model ' model ' does not exist']);
            end
            
            if isempty(deltaStr)
                m.deltaP = zeros(1,length(m.pLabels));
            else
                m.deltaP = zeros(1,length(m.pLabels));

                if ischar(deltaStr)
                    if any(strfind(deltaStr,'bias'))
                        m.deltaP(1:2) = 1;
                    end
                    
                    if any(strfind(deltaStr,'sens'))
                        m.deltaP(3:4) = 1;
                    end
                    
                    if any(strfind(deltaStr,'shape')) 
                        m.deltaP(5:end) = 1;
                    end
                    
                    if any(strfind(deltaStr,'^n'))
                        m.deltaP(5:6) = 1;
                    end
                    
                    if any(strfind(deltaStr,'^c50'))
                        m.deltaP(7:8) = 1;
                    end
                    
                    if any(strfind(deltaStr,'all'))
                        m.deltaP(1:end) = 1;
                    end
                else
                    m.deltaP = deltaStr;
                end
            end
        end

        function obj = fit(obj,behavModel,perturbationModel)
            %             f=figure;
            obj.model = obj.getModel(behavModel,perturbationModel);
            obj.lambda = 0.0001;
                       
            numSubjects = length(obj.names);
            obj.fitData.params = cell(1,numSubjects);
            obj.fitData.deltaParams = cell(1,numSubjects);
            obj.fitData.deltaCondition = cell(1,numSubjects);
            for n = 1:numSubjects
                
                idx = obj.data{n}.adapterContrast == 0.05;
                cont = obj.data{n}.stimulus(idx,:);
                resp = obj.data{n}.response(idx);
                p_nL = obj.util_optim(obj.model.ZL,obj.model.ZR,[],cont,resp,obj.model.LB,obj.model.UB);
                
                %If the base model does not fit shape params on each
                %side, make sure to copy the parameter used for each side
                p_nL(obj.model.LB==0 & obj.model.UB==0) = p_nL(find(obj.model.LB==0 & obj.model.UB==0) - 1);
                
                obj.fitData.params{n} = p_nL;
                
                %Then go through each perturbation condition and fit
                %extended models
                
                %Force UB and LB to be 0 for parameters which are NOT
                %changed by the perturbations
                LB = -inf(1,length(obj.model.pLabels));
                UB = inf(1,length(obj.model.pLabels));
                LB(~obj.model.deltaP)=0;
                UB(~obj.model.deltaP)=0;
                
                congruence = unique(obj.data{n}.adapterCongruence);
                p_L = [];
                for o = 1:length(congruence)
                    idx = obj.data{n}.adapterContrast==1 & obj.data{n}.adapterCongruence==congruence(o);
                    cont = obj.data{n}.stimulus(idx,:);
                    resp = obj.data{n}.response(idx);
                    offset = p_nL;
                    p_L(o,:) = obj.util_optim(obj.model.deltaZL,obj.model.deltaZR,offset,cont,resp,LB,UB);
                end
                
                obj.fitData.deltaParams{n} = p_L;
                obj.fitData.deltaCondition{n} = {'adapterCongruence',congruence};
                
            end

            disp('Done!'); 
            
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
                nL = obj.data{n}.adapterContrast==0.05;
                
                cont = obj.data{n}.stimulus(nL,:);
                resp = obj.data{n}.response(nL);
                
                p_nL = obj.util_optim(base.ZL,base.ZR,[],cont,resp,base.LB,base.UB);
                
                
                %If the base model does not fit shape params on each
                %side, make sure to copy the parameter used for each side
                p_nL(base.LB==0 & base.UB==0) = p_nL(find(base.LB==0 & base.UB==0) - 1);
                
                
                %Then go through each perturbation condition and fit
                %extended models
                
                cond = obj.fitData.deltaCondition{n}{2};
                
                for o = 1:length(cond)
                    L_idx = obj.data{n}.adapterContrast==1 & obj.data{n}.(obj.fitData.deltaCondition{n}{1}) == cond(o);
                    
                    cont = obj.data{n}.stimulus(L_idx,:);
                    resp = obj.data{n}.response(L_idx);
                    offset = repmat(p_nL,length(resp),1);
                    
                    loglik=[];
                    for id = 1:2
                        LB = -inf(1,length(pModel{id}.pLabels));
                        UB = inf(1,length(pModel{id}.pLabels));
                        LB(~pModel{id}.deltaP)=0;
                        UB(~pModel{id}.deltaP)=0;
                        
                        p_L = obj.util_optim(pModel{id}.deltaZL,pModel{id}.deltaZR,offset,cont,resp,LB,UB);
                        phat = obj.calculatePhat(pModel{id}.deltaZL,pModel{id}.deltaZR,offset,p_L,cont);
                        p = (resp==1).*phat(:,1) + (resp==2).*phat(:,2) + (resp==3).*phat(:,3);
                        
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
                    disp([obj.fitData.deltaCondition{n}{1} ' = ' num2str(cond(o)) '. LogLiks:[' num2str(loglik(1)) ',' num2str(loglik(2)) '] ' symbol]);
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
                idx = obj.data{n}.adapterContrast == 0.05;
                
                cont = obj.data{n}.stimulus(idx,:);
                resp = obj.data{n}.response(idx);
                
                cVals = unique(cont,'rows');
                
                numDeltaConditions = size(obj.fitData.deltaParams{n},1);
                subplot(2,numDeltaConditions+1,1);
                hold on; title('No perturbation');
                
                ph=[];
                for c = 1:size(cVals,1)
                    r = resp(cont(:,1)==cVals(c,1) & cont(:,2)==cVals(c,2));
                    [ph(c,:),pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
                    
                    dC = cVals(c,2)-cVals(c,1);
                    for ch=1:3
                        l(1)=line([1 1]*dC,pci(ch,:));
                        %                             l(2)=line([dC-0.03 dC+0.03],[1 1]*pci(ch,1));
                        %                             l(3)=line([dC-0.03 dC+0.03],[1 1]*pci(ch,2));
                        set(l,'Color',cols{ch},'Linewidth',0.5);
                    end
                end
                plot(cVals(:,2)-cVals(:,1),ph,'.','markersize',15);
                
                %Plot predictions
                testCont = [linspace(1,0,100)' zeros(100,1); zeros(100,1) linspace(0,1,100)'];
                p_hat = obj.calculatePhat(obj.model.ZL,obj.model.ZR,[],obj.fitData.params{n},testCont);
                set(gca,'ColorOrderIndex',1);
                plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
                ylim([0 1]);
                set(gcf,'color','w');
                
                
                
                %PLOT PERTURBATION DATASETS
                for delta = 1:numDeltaConditions
                    subplot(2,numDeltaConditions+1,1+delta);
                    hold on; title({[obj.fitData.deltaCondition{n}{1} ' = ' num2str(obj.fitData.deltaCondition{n}{2}(delta))],['perturbation params: ' obj.model.name{2}]});
                    
                    idx = obj.data{n}.adapterContrast == 1 & obj.data{n}.(obj.fitData.deltaCondition{n}{1}) == obj.fitData.deltaCondition{n}{2}(delta);
                    cont = obj.data{n}.stimulus(idx,:);
                    resp = obj.data{n}.response(idx);
                    cVals = unique(cont,'rows');
                    ph=[];
                    for c = 1:size(cVals,1)
                        r = resp(cont(:,1)==cVals(c,1) & cont(:,2)==cVals(c,2));
                        [ph(c,:),pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
                        
                        dC = cVals(c,2)-cVals(c,1);
                        for ch=1:3
                            l(1)=line([1 1]*dC,pci(ch,:));
                            %                             l(2)=line([dC-0.03 dC+0.03],[1 1]*pci(ch,1));
                            %                             l(3)=line([dC-0.03 dC+0.03],[1 1]*pci(ch,2));
                            set(l,'Color',cols{ch},'Linewidth',0.5);
                        end
                    end
                    plot(cVals(:,2)-cVals(:,1),ph,'.','markersize',15);
                    
                    %Plot predictions
                    testCont = [linspace(1,0,100)' zeros(100,1); zeros(100,1) linspace(0,1,100)'];
                    p_hat = obj.calculatePhat(obj.model.deltaZL,obj.model.deltaZR,obj.fitData.params{n},obj.fitData.deltaParams{n}(delta,:),testCont);
                    set(gca,'ColorOrderIndex',1);
                    plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
                    ylim([0 1]);
                    
                    
                    %PLOT PARAMETER PERTURBATIONS
                    subplot(2,numDeltaConditions+1,1+delta+numDeltaConditions+1);
                    h=bar(obj.fitData.deltaParams{n}(delta,:));
                    h.EdgeColor='none';
                    set(gca,'XTickLabel',obj.model.pLabels,'XTickLabelRotation',45);
                    
                end
            end
            
            
            
            
        end
 
        
    end
    
    methods (Access=private)
     
        function [p,exitflag] = util_optim(obj,ZL,ZR,offset,cont,resp,LB,UB)
%             lambda = 0.0001;
%             lambda=0;
            objective = @(PARAMETERS) (-obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) + obj.lambda*sum(abs(PARAMETERS).^2) );
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000,'Display','off');
            [p,~,exitflag] = fmincon(objective,zeros(1,length(LB)), [], [], [], [], LB, UB, [], options);
            %
            %             tab = tabulate(resp);tab = tab(:,3)/100;
            %             ll_0=sum(tab.*log2(tab));
            %
            %             obj.fitData.nonLaserPseudoR2{n}(s) = 1 - (-ll)/ll_0; %McFadden pseudo r^2
            %
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
        
        function logLik = calculateLogLik(obj,testParams,ZL,ZR,offset,inputs,responses)
            phat = obj.calculatePhat(ZL,ZR,offset,testParams,inputs);
            p = (responses==1).*phat(:,1) + (responses==2).*phat(:,2) + (responses==3).*phat(:,3);
%             p(p<0)=0; %Correct for weird numerical error at extreme ends of probabilities;
            %             if any(p<0)
            %                 keyboard;
            %             end
            
            p(p<0.0001) = 0.0001;
            logLik = nanmean(log2(p));
            
%             if isinf(logLik)
%                 keyboard;
%             end
            
            %             disp(logLik);
        end
        
        function phat = calculatePhat(obj,ZL,ZR,offset,testParams,inputs)
            zl = ZL(offset,testParams,inputs);
            zr = ZR(offset,testParams,inputs);
            pL = exp(zl)./(1+exp(zl)+exp(zr));
            pR = exp(zr)./(1+exp(zl)+exp(zr));
            pNG = 1 - pL - pR;

            
            phat = [pL pR pNG];
            
        end
        
        
    end
end