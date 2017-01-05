classdef nestedLogit
    properties (Access=public)
        expRef;
        modelString;
        parameterLabels;
        parameterFits;
        parameterBounds;
        parameterStart;
        Zinput;
        ZGO;
        ZLR;
        regularise;
        data;
        p_hat;
        ContrastDimensions;
    end
    
    properties (Access=public)
        guess_bpt;
    end
    
    methods
        function obj = nestedLogit(inputData)
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
                error('nestedLogit:constructorFail', 'Must pass either an expRef or data struct to GLM constructor');
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
                

                case 'C50'
                    obj.parameterLabels = {'Offset_G','Scale_G','Offset_{LR}','Scale_{LR}','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf +inf +inf];
                    
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    
                    obj.ZGO = @(P,in)( P(1) + P(2)*max( (in.^P(5))./(P(6)^P(5) + in.^P(5)) ,[],2) );
                    obj.ZLR = @(P,in)( P(3) + P(4)*diff( (in.^P(5))./(P(6)^P(5) + in.^P(5)) ,[],2) );
                    
                case 'C50-biasAsContrast'
                    obj.parameterLabels = {'Offset_G','Scale_G','Offset_{LR}','Scale_{LR}','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf +inf +inf];
                    
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    
                    obj.ZGO = @(P,in)( P(2)*(max( (in.^P(5))./(P(6)^P(5) + in.^P(5)) ,[],2) + P(1)) );
                    obj.ZLR = @(P,in)( P(4)*(diff( (in.^P(5))./(P(6)^P(5) + in.^P(5)) ,[],2) + P(3)) );
                    
                otherwise
                    error('Model does not exist');
                    
            end
           
            if isempty(obj.parameterStart)
                obj.parameterStart = zeros(1,length(obj.parameterLabels));
            end

        end
        
        function obj = fit(obj)    
            %Non crossvalidated fitting
            
            if isempty(obj.ZGO)
                error('Please set a model first using method setModel(...)');
            end
                 
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
            
            responses = obj.data.response;
            
            if isempty(obj.regularise)
                objective = @(b) (obj.calculateLogLik(b, obj.Zinput(obj.data), responses));
            else
                objective = @(b) (obj.calculateLogLik(b, obj.Zinput(obj.data), responses) + obj.regularise(b));
            end
            
            [obj.parameterFits,~,exitflag] = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
            if ~any(exitflag == [1,2])
                obj.parameterFits = nan(1,length(obj.parameterLabels));
            end
        end

        function [obj,varargout] = fitCV(obj,varargin)
            %Crossvalidated fitting
            
            if isempty(obj.ZGO)
                error('Please set a model first using method setModel(...)');
            end

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
                
                    [obj.parameterFits(f,:),~,exitflag] = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
                    if ~any(exitflag == [1,2])
                        obj.parameterFits(f,:) = nan(1,length(obj.parameterLabels));
                    end

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
                    
                    if length(uniqueC1D)<50
                        errorbar(repmat(uniqueC1D,1,r),prop(:,1:r),err(:,[1:r]),err(:,[4:(3+r)]),'.','MarkerSize',20);
                    else
                        plot(repmat(uniqueC1D,1,r),prop(:,1:r),'o');
                    end
                    xlabel('contrast');
                    ylabel('P( choice | contrast)');
                    
                    set(gca,'box','off');
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
%                         p_plot = prop(:,:,i);
%                         p_plot = [p_plot; nan(1,length(uniqueCR))];
%                         p_plot = [p_plot, nan(length(uniqueCL)+1,1)];
%                         pcolor([uniqueCR; 1],[uniqueCL; 1],p_plot); caxis([0 1]); shading('flat');
%                         set(gca,'box','off');
%                         set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);
                        imagesc(prop(:,:,i)); caxis([0 1]);
                        set(gca,'xtick',1:length(uniqueCR),'ytick',1:length(uniqueCL),'xticklabels',uniqueCR,'yticklabels',uniqueCL);
                        set(gca,'YDir','normal','box','off');
%                         xlabel('Contrast right');
%                         ylabel('Contrast left');
                        title(titles{i});
                        axis square;
%                         set(gca,'XTick','','YTick',0:0.1:0.5);
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
                        
%                         evalCL = linspace(0,max(obj.data.contrast_cond(:,1)),100);
%                         evalCR = linspace(0,max(obj.data.contrast_cond(:,1)),100);
                        evalCL = linspace(0,1,100);
                        evalCR = linspace(0,1,100);
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
                            imagesc(evalCR,evalCL,prop(:,:,i)); caxis([0 1]);
                            set(gca,'YDir','normal','box','off','YTick',0:0.2:1);
                            
                            
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
                        
                        figure('color','w');
                        cont = obj.data.contrast_cond;
                        resp = obj.data.response;
                        cVals = unique(cont(:));
%                         numPedestals = length(cVals)-1;
                        numPedestals = length(cVals);
                        cols = [0 0.4470 0.7410;
                            0.8500 0.3250 0.0980;
                            0.4940    0.1840    0.5560];
                        for ped = 1:numPedestals
                            subplot(numPedestals,1,ped); hold on;
                            
                            set(gca,'colororder',cols);
                            
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
                                    set(l,'Color',cols(ch,:),'Linewidth',0.5);
                                    %                             l.Color = cols{ch};
                                    %                             l.LineWidth=1;
                                end
                            end
                            set(gca,'ColorOrderIndex',1);
                            
                            plot(uC,ph,'.','markersize',15);
                            
                            %Plot predictions
                            testCont = [linspace(max(abs(uC))+0.1,0,100)' zeros(100,1); zeros(100,1) linspace(0,max(abs(uC))+0.1,100)'] + cVals(ped);
                            p_hat = obj.calculatePhat(obj.parameterFits,testCont);
                            set(gca,'ColorOrderIndex',1);
                            h=plot(diff(testCont,[],2),p_hat,'linewidth',1);
                            if ped == 1
                                legend(h,{'L','R','NG'}); legend boxoff;
                            end
                            xlim([-1 1]*(max(cVals)+0.1)); ylim([0 1]);
                            
                            
                            ylabel(['ped: ' num2str(cVals(ped))]);
                            
                            
                            if ped ~= numPedestals
                                set(gca,'xtick','','xcolor','w');
                            end
                            
                        end
                        
                        
                end
            else
                error('Model not fitted (non-crossvalidated) yet');
            end
        end
        

        function phat = calculatePhat(obj,testParams,inputs)
            if isempty(obj.ZGO)
                error('Please set a model first using method setModel(...)');
            end
            
            zgo = obj.ZGO(testParams,inputs);
            zlr = obj.ZLR(testParams,inputs);
            
            pNG = 1./(1+exp(zgo)); pG = 1-pNG;
            pLg = exp(zlr)./(1+exp(zlr));
            pRg = 1./(1+exp(zlr));
            
            pL = pLg .* pG;
            pR = pRg .* pG;
            
            phat = [pL pR pNG];

        end

        function plotPspace_mnrvsnested(obj)
            %use mnrfit to fit rudimentary models MNR and NESTED and
            %compare directly
            cont = obj.data.contrast_cond;
            resp = obj.data.response;
            resp_nes = resp; resp_nes(resp_nes==3)=0; resp_nes=resp_nes+1;

            B_mnr = mnrfit(cont,resp);
            B_nes = mnrfit(cont,resp_nes,'model','hierarchical');
            
            cVals = linspace(0,0.54,1000);
            [cr,cl]=meshgrid(cVals);
            cont = [cl(:) cr(:)];
            
            P_mnr = mnrval(B_mnr,cont);
            P_nes = mnrval(B_nes,cont,'model','hierarchical');
            P_nes = [P_nes(:,2:3) P_nes(:,1)];
            
            %Plot
            figure('color','w');
            labels = {'pL','pR','pNG'};
%             cVals=log(cVals);
            for r = 1:3
                subplot(1,3,r);
                [~,ax]=contour(cVals,cVals,reshape(P_mnr(:,r),size(cl)));
                ax.LineWidth=1; hold on;
                [~,ax]=contour(cVals,cVals,reshape(P_nes(:,r),size(cl)));
                ax.LineWidth=2;
                ax.LineStyle=':';
                
                set(gca,'ydir','normal','box','off'); 
                xlabel('CR'); ylabel('CL'); title(labels{r}); axis square; caxis([0 1]);
%                 
% %                 subplot(2,3,r+3);
%                 imagesc(cVals,cVals,reshape(P_nes(:,r),size(cl)));
%                 set(gca,'ydir','normal'); 
%                 xlabel('CR'); ylabel('CL'); title(labels{r});axis square;  caxis([0 1]);
            end
%             keyboard;
        end
        
  
    end
    

    
    methods (Access= {?GLM})        
        function logLik = calculateLogLik(obj,testParams, inputs, responses)
            phat = obj.calculatePhat(testParams, inputs);
            logLik = -sum(log2( phat(sub2ind(size(phat), [1:length(responses)]', responses)) ));
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