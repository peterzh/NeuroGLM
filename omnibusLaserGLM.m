classdef omnibusLaserGLM
    properties
        names = {'Spemann','Murphy','Morgan','Whipple'};
        expRefs;
        inactivationCoords;
        data;
        laserParameters;
    end
    
    properties (Access=public)
        cfn_parameters;
        moreData = 0;
        guess_bpt = [];

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
            
            cfn = @(c)(c.^obj.cfn_parameters{subj}(1))./(c.^obj.cfn_parameters{subj}(1) + obj.cfn_parameters{subj}(2).^obj.cfn_parameters{subj}(1));
            switch(model)
                case 'full_noLaser'
                    ZL = @(off,p,c)(p(1) + p(3).*cfn(c(:,1)) + p(5).*cfn(c(:,2)) );
                    ZR = @(off,p,c)(p(2) + p(4).*cfn(c(:,2)) + p(6).*cfn(c(:,1)) );
                    numP = 6;
                case 'sub_noLaser'
                    ZL = @(off,p,c)(p(1) + p(3).*cfn(c(:,1)));
                    ZR = @(off,p,c)(p(2) + p(4).*cfn(c(:,2)));
                    numP = 4;
                case 'sub_+bias'
                    ZL = @(off,p,c)(off + p(1) );
                    ZR = @(off,p,c)(off + p(2) );
                    numP = 2;
                case 'sub_+sens'
                    ZL = @(off,p,c)(off + p(1).*cfn(c(:,1)) );
                    ZR = @(off,p,c)(off + p(2).*cfn(c(:,2)) );
                    numP = 2;
                case 'sub_+bias+sens'
                    ZL = @(off,p,c)(off + p(1) + p(3).*cfn(c(:,1)) );
                    ZR = @(off,p,c)(off + p(2) + p(4).*cfn(c(:,2)) );
                    numP = 4;
            end
        end
            
        
        function obj = fit(obj,varargin)
            nonLaserModel = 'sub_noLaser';
            laserModel = 'sub_+bias+sens';
            numSubjects = length(obj.names);
            obj.laserParameters = cell(1,numSubjects);
            for n = 1:numSubjects

                options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
                
                %First fit non-laser portion of the data
                [ZL_nL,ZR_nL,numP] = obj.getModel(nonLaserModel,n);
                nL = obj.data{n}.laserIdx==0;

                cont = obj.data{n}.contrast_cond(nL,:); 
                resp = obj.data{n}.response(nL); 
                
                offset = zeros(length(resp),2);
                objective = @(PARAMETERS) (obj.calculateLogLik(PARAMETERS,ZL_nL,ZR_nL,offset,cont,resp));
                
                options = optiset('display','final','solver','NOMAD');
                Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                [p_nL,~,exitflag,~] = Opt.solve;
                if exitflag < 0
                    error('Bad fit');
                end
                
                %Then go through each inactivation location and fit
                %extended models
                p_L = [];
                for site = 1:size(obj.inactivationCoords{n},1)
                    L_idx = obj.data{n}.laserIdx==site;
                    
                    cont = obj.data{n}.contrast_cond(L_idx,:);
                    resp = obj.data{n}.response(L_idx);
                    
                    offset = [ZL_nL([],p_nL,cont) ZR_nL([],p_nL,cont)];
                    
                    [ZL,ZR,numP] = obj.getModel(laserModel,n);
                    objective = @(PARAMETERS) (obj.calculateLogLik(PARAMETERS,ZL,ZR,offset,cont,resp));
                    Opt = opti('fun',objective,'x0',zeros(1,numP),'options',options);
                    p_L(site,:) = Opt.solve';
                end                
                obj.laserParameters{n} = p_L;
            end
            disp('Done!');
        end
        
        function plotModel(obj)
            numSubjects = length(obj.names);
            figure;
            
            param_labels = {'b_L','b_R','sens_L','sens_R'};
            for n = 1:numSubjects
                p_L = obj.laserParameters{n};
                for i = 1:size(p_L,2)
                    subplot(numSubjects,4,4*n-4 + i);
                    scatter(obj.inactivationCoords{n}(:,2),obj.inactivationCoords{n}(:,1),200,p_L(:,i),'s','filled'); axis equal; caxis([-1 1]*2);
                    set(gca,'box','off','ytick','','xtick','');
                    ylabel(obj.names{n});
                    
                    if n==1
                        title(param_labels{i});
                    end
                end
                
                cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                colormap(cmap);
            end
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
        
        
    end
    
    methods (Access=private)
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
            logLik = -sum(log(p));
        end
        
        function phat = calculatePhat(obj,ZL,ZR,offset,testParams,inputs)            
            zl = ZL(offset(:,1),testParams,inputs);
            zr = ZR(offset(:,2),testParams,inputs);
            pL = exp(zl)./(1+exp(zl)+exp(zr));
            pR = exp(zr)./(1+exp(zl)+exp(zr));
            pNG = 1 - pL - pR;
            
            phat = [pL pR pNG];
            
        end
        
    end
end