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
    
    properties (Access=private)
        cfn_parameters;
        moreData = 0;
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
                case 'sub_+bias'
                    ZL = @(offL,p,c)(offL + p(1) );
                    ZR = @(offR,p,c)(offR + p(2) );
                    numP = 2;
                case 'sub_+sens'
                    ZL = @(offL,p,c)(offL + p(1).*cfn(c(:,1)) );
                    ZR = @(offR,p,c)(offR + p(2).*cfn(c(:,2)) );
                    numP = 2;
                case 'sub_+bias+sens'
                    ZL = @(offL,p,c)(offL + p(1) + p(3).*cfn(c(:,1)) );
                    ZR = @(offR,p,c)(offR + p(2) + p(4).*cfn(c(:,2)) );
                    numP = 4;
                case 'sub_+bias+sens_matchedUnits'
                    ZL = @(offL,p,c)(offL + p(1)*p(3) + p(3).*cfn(c(:,1)) );
                    ZR = @(offR,p,c)(offR + p(2)*p(4) + p(4).*cfn(c(:,2)) );
                    numP = 4;
            end
        end
                   
        function obj = fit(obj)
            f=figure;
            obj.fitData.laserModel = 'sub_+bias+sens';
            obj.fitData.lambda = 'lambda_min';

            numSubjects = length(obj.names);
%             numSubjects = 1;
            obj.fitData.params = cell(1,numSubjects);
            for n = 1:numSubjects
               
                %First fit non-laser portion of the data
%                 [ZL_nL,ZR_nL,numP] = obj.getModel(nonLaserModel,n);
                nL = obj.data{n}.laserIdx==0;

                cont = obj.data{n}.contrast_cond(nL,:); 
                resp = obj.data{n}.response(nL); 
                X = obj.getNonLaserDesignMatrix(n,obj.cfn(cont,n),obj.data{n}.sessionID(nL));
                opts = []; opts.intr=0; %exclude bias
                fit = cvglmnet(X,resp,'multinomial',glmnetSet(opts));
                
                %Then go through each inactivation location and fit
                %extended models
                p_L = [];
                for site = 1:size(obj.inactivationCoords{n},1)
                    L_idx = obj.data{n}.laserIdx==site;
                    
                    cont = obj.data{n}.contrast_cond(L_idx,:); 
                    X = obj.getNonLaserDesignMatrix(n,obj.cfn(cont,n),obj.data{n}.sessionID(L_idx));
                    resp = obj.data{n}.response(L_idx);
                    
%                     offset = [ZL_nL([],p_nL,cont) ZR_nL([],p_nL,cont)];
                    %Calculate offset from non-laser model fit
                    offset = cvglmnetPredict(fit,X,obj.fitData.lambda,'link'); offset=[offset(:,1)-offset(:,3) offset(:,2)-offset(:,3)];
                    
                    %Use that offset to manually fit a model at each
                    %inactivation site
                    [ZL,ZR,numP] = obj.getModel(obj.fitData.laserModel,n);
                    objective = @(PARAMETERS) (obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp) + 0.01*sum(PARAMETERS.^2));
%                     keyboard;
                    LB = ones(1,numP)*-10;
                    UB = ones(1,numP)*10;
                    
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'bounds',LB,UB,'options',options);
                            p_L(site,:) = Opt.solve';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                            p_L(site,:) = fmincon(objective,zeros(1,numP), [], [], [], [], LB, UB, [], options);
                    end 
                    
                    obj.fitData.siteX{site,n} = X;
                    obj.fitData.siteFit{site,n} = fit;
                end
                
%                 keyboard;
                obj.fitData.params{n} = p_L;
            end

            
            disp('Done!'); close(f);
        end
        
        function obj = fitNull(obj)
            laserModel = 'sub_+bias+sens';
            lambda = 0;
            numSubjects = length(obj.names);
%             numSubjects=1;
            numIter = 200;
            
            obj.fitData.paramsNull = cell(1,numSubjects);
            for n = 1:numSubjects
                %Get all non-laser trials
                nL = obj.data{n}.laserIdx==0;
                tab = tabulate(obj.data{n}.laserIdx);
                proportion = mean(tab(2:end,3)/100);
                
                cont = obj.data{n}.contrast_cond(nL,:); 
                resp = obj.data{n}.response(nL); 
                X = obj.getNonLaserDesignMatrix(n,obj.cfn(cont,n),obj.data{n}.sessionID(nL));
                
                for iter = 1:numIter
                    disp(iter);
                    %In each interation, assign a proportion of them to be
                    %a laser site
                    q_idx = randsample([1:size(X,1)]',round(proportion*size(X,1)));
                    
                    %Of the remaining nonLaser trials, fit the nonLaser
                    %model
                    XnL = X; XnL(q_idx,:) = [];
                    RnL = resp; RnL(q_idx) = [];
                    opts = []; opts.intr=0; %exclude bias
                    fit = cvglmnet(XnL,RnL,'multinomial',glmnetSet(opts));
                    
                    %Now fit the laser site using the offset defined from
                    %the non-laser component
                    XL = X(q_idx,:);
                    RL = resp(q_idx);
                    contL = cont(q_idx,:);
                    
                    offset = cvglmnetPredict(fit,XL,lambda,'link'); 
                    offset=[offset(:,1)-offset(:,3) offset(:,2)-offset(:,3)];
                    
                    [ZL,ZR,numP] = obj.getModel(laserModel,n);
                    objective = @(PARAMETERS) (obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,contL,RL) + 0.01*sum(PARAMETERS.^2));
                    
                    switch(obj.fitMethod)
                        case 'opti'
                            options = optiset('display','final','solver','NOMAD');
                            Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                            obj.fitData.paramsNull{n}(iter,:) = Opt.solve';
                            
                        case 'fmincon'
                            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                            obj.fitData.paramsNull{n}(iter,:) = fmincon(objective,zeros(1,numP), [], [], [], [], [], [], [], options);
                    end
                end
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
        
        function obj = plotModel(obj)
            numSubjects = length(obj.names);
            figure;
            
            param_labels = {'b_L','b_R','sens_L','sens_R'};
            for n = 1:numSubjects
                p_L = obj.fitData.params{n};
                for i = 1:size(p_L,2)
                    subplot(numSubjects,4,4*n-4 + i);
                    scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),200,p_L(:,i),'s','filled'); axis equal; 
                    
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
            
            figure;
            param_labels = {'b_L - b_R','b_L + b_R'};
            for n = 1:numSubjects
                p_L = obj.fitData.params{n};
                plotVal = {p_L(:,1)-p_L(:,2), p_L(:,1)+p_L(:,2)};
                for i = 1:2
                    subplot(numSubjects,4,2*n-2 + i);
                    scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),200,plotVal{i},'s','filled'); axis equal; 
                    
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
            switch(what)
                case 'performance'
                    figure;
                    kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
                    clear PD;
                    for s = 1:numSubjects
                        subplot(1,numSubjects,s);
                        imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
                        set(imX,'alphadata',0.7);
                        hold on;
                        
                        PF = 100*pivottable(obj.data{s}.laserIdx,[],obj.data{s}.feedbackType==1,'mean');
                        
                        % NULL TEST
                        tab = tabulate(obj.data{s}.laserIdx);
                        proportion = mean(tab(2:end,3)/100);
                        NUM_SHUFFLES = 1000;
                        nullP = nan(NUM_SHUFFLES,1);
                        nonL_idx = obj.data{s}.laserIdx==0;
                        
                        for shuff = 1:NUM_SHUFFLES
                            disp(shuff);
                            L_idx = obj.data{s}.laserIdx;
                            nL = randsample(find(nonL_idx),round(proportion*sum(nonL_idx)));
                            L_idx(nL) = max(obj.data{s}.laserIdx)+1;
                            
                            shuf_p = 100*pivottable(L_idx,[],obj.data{s}.feedbackType==1,'mean');
                            nullP(shuff) = shuf_p(end);
                        end
                        
                        nullP_ci95 = quantile(nullP,[0.025 0.975]);
                        nullP_ci99 = quantile(nullP,[0.005 0.995]);
                        
                        
                        ap = obj.inactivationCoords{s}(:,1);
                        ml = obj.inactivationCoords{s}(:,2);
                        out = PF(2:end);
                        if bilateral_mirror == 1
                            ap = [ap;ap];
                            ml = [ml;-ml];
                            out = [out;out];
                        end
                        
                        %sig testing
                        sig95 = double((out < nullP_ci95(1)) | (nullP_ci95(2) < out));
                        sig99 = double((out < nullP_ci99(1)) | (nullP_ci99(2) < out));
                        sig = sum([sig95 sig99],2);
                        sig(sig==0)=0.2;
                        
                        s1=scatter(ml,ap,sig*80,out,'s','filled'); axis equal; colorbar;
                        %         hold on;
                        s1.MarkerEdgeColor=[1 1 1]*0.9;
                        
                        %     s1=scatter(ml,ap,150,out,'s','filled'); axis equal; colorbar;
                        ylim([-6 4]); %s1.MarkerEdgeColor=[1 1 1]*0.8;
                        cax = PF(1) + [-1 1]*max(abs(PF(2:end)-PF(1)));
                        try
                            PD(:,s)=PF(2:end)-PF(1);
                        catch
                        end
                        %     cax(cax>1)=1;
                        caxis(cax);
                        hold off;
                        title(obj.names{s})
                        xlim([-5 5]);
                        
                        cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                            linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                        colormap(flipud(cmap));
                        
                        set(gca,'xtick','','ytick','');
                        set(gca,'box','off','XColor','w','YColor','w');
                        
                    end
                    
                    cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
                        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                    colormap(flipud(cmap));
                    
                    set(gca,'xtick','','ytick','');
                    set(gca,'box','off','XColor','w','YColor','w');
                    set(gcf,'Color','w');
                    
%                     figure;
                    % imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal');
                    % set(imX,'alphadata',0.7);
                    % hold on;
                    
%                     ap = l(s).inactivationSite(:,1);
%                     ml = l(s).inactivationSite(:,2);
%                     out = mean(PD,2);
%                     if bilateral_mirror == 1
%                         ap = [ap;ap];
%                         ml = [ml;-ml];
%                         out = [out;out];
%                     end
%                     
%                     s1=scatter(ml,ap,150,out,'s','filled'); axis equal; colorbar;
%                     ylim([-6 4]); %s1.MarkerEdgeColor=[1 1 1]*0.8;
%                     cax = [-1 1]*max(abs(mean(PD,2)));
%                     caxis(cax);
%                     hold off;
%                     xlim([-5 5]);
%                     cmap = [ linspace(0,1,100)' ones(100,1) linspace(0,1,100)';
%                         linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
%                     colormap(flipud(cmap));
%                     set(gca,'xtick','','ytick','');
%                     set(gca,'box','off','XColor','w','YColor','w');
%                     
%                     set(gcf,'Color','w');
                case 'RT'
                case 'choices'
            end
        end
        
        function plotLogLik(obj)
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
                    PLOTRANGES = [-1 1 -3 3]*2;
                    SAMPLING = 100;
%                     PLOTRANGES = [-1 1]*max(abs(p_L(:)));
                    ax=axes; hold on;
                    c1=ezcontourf(fcn_L,PLOTRANGES,SAMPLING);   
                    
%                     bL_vals = linspace(-5,5,SAMPLING);
%                     sL_vals = linspace(-15,15,SAMPLING);
%                     
%                     img = nan(length(sL_vals),length(bL_vals));
%                     for bl = 1:length(bL_vals)
%                         for sl = 1:length(sL_vals)
%                             img(sl,bl) = fcn_L(bL_vals(bl),sL_vals(sl));
%                         end
%                     end
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

            imagesc(X); drawnow;
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
            
            logLik = -sum(log(p));
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