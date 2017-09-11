classdef omnibusLaserGLM
    properties
        expt;
        names;
        expRefs;
        inactivationCoords;
        data;
        data_shuffle;
        model;
        lambda;
        fitData;
        cvData;
        guess_bpt = [];
    end
    
    properties (Access=public)
        %         cfn_parameters;
        significance_flag = 0;
        bilateral_flag = 0;
        moreData = 0;
        nestedModel_flag = 0;
    end
    
    methods
        function obj = omnibusLaserGLM(expt,nameList)
            %             Constructor which takes in a list of subject names (cell array),
            %             and loads the data
            
            disp('Collating expRefs...');
            obj.names = nameList;
            obj.expt = expt;
            numSubjects = length(obj.names);
            obj.expRefs = cell(1,numSubjects);
            for n = 1:numSubjects
                eRefs = obj.getExpRefs(obj.names{n},expt);
                obj.expRefs{n} = eRefs([eRefs{:,2}]'==1,1);
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
                    d = loadData(obj.expRefs{n}{session});
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
            
            if contains(expt,'bilateral')
                obj.bilateral_flag = 1;
                warning('SET BILATERAL FLAG');
            end
            
            %match laserIdx over all subjects
            %create laserIdx variable
            allD = [obj.data{:}];
            allLas = cat(1,allD.laser);
            [sites,~,ic] = unique(allLas,'rows');
            nanIdx = find(isnan(sites(:,1)),1,'first');
            inactivationSite = sites(1:nanIdx-1,:);
            ic(ic>=nanIdx)=0;
            
            N=arrayfun(@(s)length(s.response),allD);
            for n = 1:numSubjects
                obj.data{n}.laserIdx = ic(1:N(1));
                ic(1:N(1)) = [];
                N(1)=[];
                
                obj.data{n}.areaIdx = obj.identifyArea(obj.data{n}.laser);
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
            
            %             numSubjects = numSubjects + 1;
            
            %             %now fit c50 and n parameters for each subject's nonlaser data
            %             obj.cfn_parameters = cell(1,numSubjects);
            %             for n = 1:numSubjects
            %                 g=GLM(getrow(obj.data{n},obj.data{n}.laserIdx==0)).setModel('C50-subset').fit;
            %                 obj.cfn_parameters{n} = g.parameterFits(5:6);
            %             end
            %             obj.plotExptSummary;
            
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
        
        function m = getModel2(~,model,deltaStr)
            %More efficient model definitions
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
                    m.LB = [-1 -1 0 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+1 +1 +inf +inf 1 0]; %only used when doing non-laser fit
                    m.pLabels = {'bL','bR','sL','sR','nL','nR'};
                    
                    m.deltaZL = @(p_offset,p,c)( 10*(p_offset(:,3)*2^p(3)).*(m.cfn(c(:,1),p_offset(:,5)*2^p(5)) + p_offset(:,1) + p(1)) );
                    m.deltaZR = @(p_offset,p,c)( 10*(p_offset(:,4)*2^p(4)).*(m.cfn(c(:,2),p_offset(:,6)*2^p(6)) + p_offset(:,2) + p(2)) );
                    
                case 'BO_NC50'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p_offset,p,c)( p(1) + 10*p(3).*m.cfn(c(:,1),p(5),p(7)) );
                    m.ZR = @(p_offset,p,c)( p(2) + 10*p(4).*m.cfn(c(:,2),p(5),p(7)) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','sL','sR','nL','nR','c50L','c50R'};
                    
                    m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,3)*2^p(3)).*m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) );
                    m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + 10*(p_offset(:,4)*2^p(4)).*m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8))  );
                    
                case 'BO_NC50_DISTORTION'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p_offset,p,c)( p(1) + 10*p(3).*( m.cfn(c(:,1),p(5),p(7)) - m.cfn(c(:,2),p(5),p(7)).*(m.cfn(c(:,1),p(5),p(7)) - p(9))*1/3 ) );
                    m.ZR = @(p_offset,p,c)( p(2) + 10*p(4).*( m.cfn(c(:,2),p(5),p(7)) - m.cfn(c(:,1),p(5),p(7)).*(m.cfn(c(:,2),p(5),p(7)) - p(10))*1/3 ) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0 1 1]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','sL','sR','nL','nR','c50L','c50R','phiL','phiR'};
                    
                    m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,3)*2.^p(3)).*( m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) - m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)).*(m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) - (p_offset(:,9)+p(9)) )*1/3 ) );
                    m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + 10*(p_offset(:,4)*2.^p(4)).*( m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)) - m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)).*(m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)) - (p_offset(:,10)+p(10)) )*1/3 ) );
                    
                case 'BO_NC50_NORMALISATION'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p_offset,p,c)( p(1) + 10*p(9).*m.cfn(c(:,1),p(5),p(7)) + 10*p(3).*m.cfn(c(:,1),p(5),p(7))./(m.cfn(c(:,1),p(5),p(7)) + m.cfn(c(:,2),p(5),p(7))) );
                    m.ZR = @(p_offset,p,c)( p(2) + 10*p(10).*m.cfn(c(:,2),p(5),p(7)) + 10*p(4).*m.cfn(c(:,2),p(5),p(7))./(m.cfn(c(:,1),p(5),p(7)) + m.cfn(c(:,2),p(5),p(7))) );
                    m.LB = [-inf -inf 0 0 0.3 0 0.001 0 0 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0 +inf +inf]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','snL','snR','nL','nR','c50L','c50R','sL','sR'};
                    
                    m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,3)*2.^p(3)).*m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7))./(m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) + m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)))  );
                    m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + 10*(p_offset(:,4)*2.^p(4)).*m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8))./(m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) + m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)))  );
                    
                    
                case 'BO_NC50_2SIDE'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p_offset,p,c)( p(1) + 10*p(3).*m.cfn(c(:,1),p(5),p(7)) + 10*p(5).*m.cfn(c(:,2),p(5),p(7)) );
                    m.ZR = @(p_offset,p,c)( p(2) + 10*p(4).*m.cfn(c(:,2),p(5),p(7)) + 10*p(6).*m.cfn(c(:,1),p(5),p(7)) );
                    m.LB = [-inf -inf 0 0 0 0 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'oL','oR','sL','sR','sL.','sR.','nL','nR','c50L','c50R'};
                    
                    m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,3)*2^p(3)).*m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) + 10*(p_offset(:,5)*2^p(5)).*m.cfn(c(:,2),p_offset(:,5),p_offset(:,7)) );
                    m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + 10*(p_offset(:,4)*2^p(4)).*m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)) + 10*(p_offset(:,6)*2^p(6)).*m.cfn(c(:,1),p_offset(:,6),p_offset(:,8)) );
                    
                case 'BC_NC50'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    %                     m.cfn = @(c,varargin)((c.^varargin{1})./(c.^varargin{1} + varargin{2}.^varargin{1}));
                    m.ZL = @(p_offset,p,c)( 10*p(3).*(m.cfn(c(:,1),p(5),p(7)) + p(1)) );
                    m.ZR = @(p_offset,p,c)( 10*p(4).*(m.cfn(c(:,2),p(5),p(7)) + p(2)) );
                    m.LB = [-1 -1 0 0 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+1 +1 +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'bL','bR','sL','sR','nL','nR','c50L','c50R'};
                    
                    m.deltaZL = @(p_offset,p,c)( 10*(p_offset(:,3)*2^p(3)).*(m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) + p_offset(:,1) + p(1)) );
                    m.deltaZR = @(p_offset,p,c)( 10*(p_offset(:,4)*2^p(4)).*(m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)) + p_offset(:,2) + p(2)) );
                    
                case 'BO_NC50_NESTED'
                    m.cfn = @(c,n,c50)((c.^n)./(c.^n + c50.^n));
                    m.ZL = @(p_offset,p,c)( p(1) + p(3).*max([m.cfn(c(:,1),p(5),p(7)) m.cfn(c(:,2),p(5),p(7))],[],2) ); %actually ZGnG
                    m.ZR = @(p_offset,p,c)( p(2) + p(4).*diff([m.cfn(c(:,1),p(5),p(7)) m.cfn(c(:,2),p(5),p(7))],[],2) ); %actually ZLR
                    m.LB = [-inf -inf -inf -inf 0.3 0 0.001 0]; %only used when doing non-laser fit
                    m.UB = [+inf +inf +inf +inf 20 0 3 0]; %only used when doing non-laser fit
                    m.pLabels = {'oGO','oLR','sGO','sLR','nL','nR','c50L','c50R'};
                    
                    m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + (p_offset(:,3)*2^p(3)).*max([m.cfn(c(:,1),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7)) m.cfn(c(:,2),p_offset(:,5)*2^p(5),p_offset(:,7)*2^p(7))],[],2) ); %actually ZGnG
                    m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + (p_offset(:,4)*2^p(4)).*diff([m.cfn(c(:,1),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8)) m.cfn(c(:,2),p_offset(:,6)*2^p(6),p_offset(:,8)*2^p(8))],[],2) ); %actually ZLR
                    
                    
                    %                 case 'BO_C50'
                    %                     N=1;
                    %                     m.cfn = @(c,c50)((c.^N)./(c.^N + c50.^N));
                    %                     m.ZL = @(p_offset,p,c)( p(1) + 10*p(3).*m.cfn(c(:,1),p(5)) );
                    %                     m.ZR = @(p_offset,p,c)( p(2) + 10*p(4).*m.cfn(c(:,2),p(5)) );
                    %                     m.LB = [-inf -inf 0 0 0.001 0]; %only used when doing non-laser fit
                    %                     m.UB = [+inf +inf +inf +inf 3 0]; %only used when doing non-laser fit
                    %                     m.pLabels = {'oL','oR','sL','sR','c50L','c50R'};
                    %
                    %                     m.deltaZL = @(p_offset,p,c)( p_offset(:,1) + p(1) + 10*(p_offset(:,3)*2^p(3)).*m.cfn(c(:,1),p_offset(:,5)*2^p(5)) );
                    %                     m.deltaZR = @(p_offset,p,c)( p_offset(:,2) + p(2) + 10*(p_offset(:,4)*2^p(4)).*m.cfn(c(:,2),p_offset(:,6)*2^p(6))  );
                    %
                    %                 case 'BC_C50'
                    %                     N=1;
                    %                     m.cfn = @(c,c50)((c.^N)./(c.^N + c50.^N));
                    %                     m.ZL = @(p_offset,p,c)( 10*p(3).*(m.cfn(c(:,1),p(5)) + p(1)) );
                    %                     m.ZR = @(p_offset,p,c)( 10*p(4).*(m.cfn(c(:,2),p(5)) + p(2)) );
                    %                     m.LB = [-inf -inf 0 0 0.001 0]; %only used when doing non-laser fit
                    %                     m.UB = [+inf +inf +inf +inf 3 0]; %only used when doing non-laser fit
                    %                     m.pLabels = {'bL','bR','sL','sR','c50L','c50R'};
                    %
                    %                     m.deltaZL = @(p_offset,p,c)( 10*(p_offset(:,3)*2^p(3)).*(m.cfn(c(:,1),p_offset(:,5)*2^p(5)) + p_offset(:,1) + p(1)) );
                    %                     m.deltaZR = @(p_offset,p,c)( 10*(p_offset(:,4)*2^p(4)).*(m.cfn(c(:,2),p_offset(:,6)*2^p(6)) + p_offset(:,2) + p(2)) );
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
                        %                         m.deltaP(end-1:end) = 1;
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
        
        function data_shuffle = shuffleData(obj,n)
            
            %Per session, contrast condition,
            %shuffle laser location identities
            
            data_shuffle = obj.data{n};
            numSessions = max(obj.data{n}.sessionID);
            
            cont = obj.data{n}.stimulus(:,1:2);
            c_0 = sum(cont,2)==0;
            c_L = cont(:,1) > 0 & cont(:,2)==0;
            c_R = cont(:,2) > 0 & cont(:,1)==0;
            c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
            cont = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
            
            for sess = 1:numSessions
                for con = 1:4
                    idx = obj.data{n}.sessionID==sess & cont == con;
                    
                    old_lasIdx = obj.data{n}.laserIdx(idx);
                    new_lasIdx = old_lasIdx(randperm(length(old_lasIdx)));
                    
                    data_shuffle.laserIdx(idx) = new_lasIdx;
                    
                    old_areaIdx = obj.data{n}.areaIdx(idx);
                    new_areaIdx = old_areaIdx(randperm(length(old_areaIdx)));
                    
                    data_shuffle.areaIdx(idx) = new_areaIdx;
                end
            end
            
            
        end
        
        function obj = fit(obj,nonLaserModel,LaserModel)
            %             f=figure;
            obj.model = obj.getModel2(nonLaserModel,LaserModel);
            obj.fitData.perSession = 0;
            obj.fitData.groupAreas = 1;
            obj.lambda = 0.0001;
            
            if strcmp(nonLaserModel(end-5:end),'NESTED')
                obj.nestedModel_flag = 1;
            else
                obj.nestedModel_flag = 0;
            end
            
            numSubjects = length(obj.names);
            obj.fitData.params = cell(1,numSubjects);
            obj.fitData.deltaParams = cell(1,numSubjects);
            for n = 1:numSubjects
                
                %First fit non-laser portion of the data
                
                %Manual optim method per-session
                nL = obj.data{n}.laserIdx==0;
                
                sessions = unique(obj.data{n}.sessionID(nL));
                if obj.fitData.perSession==1
                    p_nL = nan(length(sessions),length(obj.model.pLabels));
                    %                     obj.fitData.nonLaserPseudoR2{n} = nan(length(sessions),1);
                    for s = 1:length(sessions)
                        cont = obj.data{n}.stimulus(nL & obj.data{n}.sessionID==sessions(s),1:2);
                        resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                        
                        p_nL(s,:) = obj.util_optim(obj.model.ZL,obj.model.ZR,[],cont,resp,obj.model.LB,obj.model.UB);
                    end
                else
                    p_nL = nan(1,length(obj.model.pLabels));
                    cont = obj.data{n}.stimulus(nL,:);
                    resp = obj.data{n}.response(nL);
                    
                    p_nL(1,:) = obj.util_optim(obj.model.ZL,obj.model.ZR,[],cont,resp,obj.model.LB,obj.model.UB);
                    p_nL = repmat(p_nL,length(sessions),1);
                end
                
                
                %If the base model does not fit shape params on each
                %side, make sure to copy the parameter used for each side
                p_nL(:,obj.model.LB==0 & obj.model.UB==0) = p_nL(:,find(obj.model.LB==0 & obj.model.UB==0) - 1);
                
                obj.fitData.params{n} = p_nL;
                
                %                 if n == 7
                %                     keyboard;
                %                 end
                
                %Then go through each inactivation location and fit
                %extended models
                if obj.fitData.groupAreas==0
                    laserIdx = 1:size(obj.inactivationCoords,1);
                    siteID = obj.data{n}.laserIdx;
                else
                    laserIdx = (1:6);
                    siteID = obj.data{n}.areaIdx;
                end
                
                
                %Force UB and LB to be 0 for laser parameters which are NOT
                %changed by the laser
                LB = -inf(1,length(obj.model.pLabels));
                UB = inf(1,length(obj.model.pLabels));
                LB(~obj.model.deltaP)=0;
                UB(~obj.model.deltaP)=0;
                
                %                 p_L = nan(size(obj.inactivationCoords,1),length(UB));
                p_L = nan(max(laserIdx),length(UB));
                
                for site = laserIdx
                    %                     L_idx = obj.data{n}.laserIdx==site;
                    L_idx = siteID==site;
                    if sum(L_idx)>0
                        cont = obj.data{n}.stimulus(L_idx,:);
                        resp = obj.data{n}.response(L_idx);
                        sess = obj.data{n}.sessionID(L_idx);
                        
                        %Calculate offset from non-laser fit
                        %offset on parameters
                        offset = p_nL(sess,:);
                        
                        p_L(site,:) = obj.util_optim(obj.model.deltaZL,obj.model.deltaZR,offset,cont,resp,LB,UB);
                    end
                end
                obj.fitData.deltaParams{n} = p_L;
            end
            
            disp('Done!'); %close(f);
            
        end
        
        function deltaShuffle(obj,nonLaserModel,LaserModel)
            m = obj.getModel2(nonLaserModel,LaserModel);
            n = 6;
            if obj.fitData.groupAreas == 0
                error('need to run grouped area fit');
            end
            
            %If the laser has no true effect, then laser trials are
            %indistinguishable from nonlaser trials
            
            %Go through each location, and shuffle with nonlaser trials,
            %and do a full fit process
            areaList = [1 2 5 6];
            numIter = 100;
            obj.lambda = 0.0001;
            
            for area = 1:length(areaList)
                D = getrow(obj.data{n},obj.data{n}.areaIdx==999 | obj.data{n}.areaIdx==areaList(area));
                
                p_L = nan(numIter,length(m.UB));
                for iter = 1:numIter
                    disp([num2str(iter) '/' num2str(numIter)]);
                    %Shuffle laser identity
                    D.areaIdx = D.areaIdx(randperm(length(D.areaIdx)));
                    
                    %Get 'nonlaser' trials, fit model
                    nL = D.areaIdx==999;
                    sessions = unique(D.sessionID(nL));
                    p_nL = nan(length(sessions),length(m.pLabels));
                    for s = 1:length(sessions)
                        cont = D.stimulus(nL & D.sessionID==sessions(s),1:2);
                        resp = D.response(nL & D.sessionID==sessions(s));
                        
                        p_nL(s,:) = obj.util_optim(m.ZL,m.ZR,[],cont,resp,m.LB,m.UB);
                    end
                    p_nL(:,m.LB==0 & m.UB==0) = p_nL(:,find(m.LB==0 & m.UB==0) - 1);
                    
                    
                    %Get 'laser' trials, fit model
                    L_idx = D.areaIdx~=999;
                    
                    LB = -inf(1,length(m.pLabels));
                    UB = inf(1,length(m.pLabels));
                    LB(~m.deltaP)=0;
                    UB(~m.deltaP)=0;
                    
                    cont = D.stimulus(L_idx,:);
                    resp = D.response(L_idx);
                    sess = D.sessionID(L_idx);
                    offset = p_nL(sess,:);
                    
                    p_L(iter,:) = obj.util_optim(m.deltaZL,m.deltaZR,offset,cont,resp,LB,UB);
                end
                
                keyboard;
            end
            
        end
        
        function likelihoodratiotest(obj,nonLaserModel,lasermodel_unrestricted,lasermodel_restricted)
            obj.lambda = 0;
            
            numSubjects = length(obj.names);
            params = cell(1,numSubjects);
            deltaParams = cell(1,numSubjects);
            
            nonlasermodel = obj.getModel2(nonLaserModel,[]);
            lasModel{1} = obj.getModel2(nonLaserModel,lasermodel_unrestricted);
            lasModel{2} = obj.getModel2(nonLaserModel,lasermodel_restricted);
            
            areaLabels = {'Lv1','Rv1','Ls1','Rs1','Lm2','Rm2'};
            
            for n = 1:numSubjects
                
                %                 if n == numSubjects
                %                     keyboard;
                %                 end
                %First fit non-laser portion of the data
                
                %Manual optim method per-session
                nL = obj.data{n}.laserIdx==0;
                
                sessions = unique(obj.data{n}.sessionID(nL));
                p_nL = nan(length(sessions),length(nonlasermodel.pLabels));
                for s = 1:length(sessions)
                    cont = obj.data{n}.stimulus(nL & obj.data{n}.sessionID==sessions(s),1:2);
                    resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                    
                    p_nL(s,:) = obj.util_optim(nonlasermodel.ZL,nonlasermodel.ZR,[],cont,resp,nonlasermodel.LB,nonlasermodel.UB);
                end
                
                %If the base model does not fit shape params on each
                %side, make sure to copy the parameter used for each side
                p_nL(:,nonlasermodel.LB==0 & nonlasermodel.UB==0) = p_nL(:,find(nonlasermodel.LB==0 & nonlasermodel.UB==0) - 1);
                
                
                %Then go through each inactivation location and fit
                %extended models
                laserIdx = (1:6);
                siteID = obj.data{n}.areaIdx;
                
                for site = laserIdx
                    L_idx = siteID==site;
                    if sum(L_idx)>0
                        cont = obj.data{n}.stimulus(L_idx,1:2);
                        resp = obj.data{n}.response(L_idx);
                        sess = obj.data{n}.sessionID(L_idx);
                        
                        %Calculate offset from non-laser fit
                        %offset on parameters
                        offset = p_nL(sess,:);
                        
                        loglik=[];
                        for id = 1:2
                            LB = -inf(1,length(lasModel{id}.pLabels));
                            UB = inf(1,length(lasModel{id}.pLabels));
                            LB(~lasModel{id}.deltaP)=0;
                            UB(~lasModel{id}.deltaP)=0;
                            
                            p_L = obj.util_optim(lasModel{id}.deltaZL,lasModel{id}.deltaZR,offset,cont,resp,LB,UB);
                            phat = obj.calculatePhat(lasModel{id}.deltaZL,lasModel{id}.deltaZR,offset,p_L,cont);
                            p = (resp==1).*phat(:,1) + (resp==2).*phat(:,2) + (resp==3).*phat(:,3);
                            
                            loglik(id) = sum(log(p));
                            
                        end
                        
                        dof = sum(lasModel{1}.deltaP) - sum(lasModel{2}.deltaP);
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
                        disp([obj.names{n} ' ' areaLabels{site} ' ' symbol]);
                        
                        
                    end
                    
                end
            end
            
        end
        
        function obj = simulateDataset(obj,profile)
            %Inserts an extra dataset with simulated nonlaser and laser
            %model behaviour
            
            if ~any(strcmp(obj.names,'<<SIMULATION>>')==1)
                simID = length(obj.names)+1;
                
                obj.names{simID}  = '<<SIMULATION>>';
                obj.data{simID} = struct;
                obj.guess_bpt(simID) = nan;
            else
                simID = length(obj.names);
            end
            
            numTrials = round(length(obj.data{simID-1}.response)/5);
            sampleIdx = randsample(length(obj.data{simID-1}.response),numTrials);
            
            D = structfun(@(field)(field(sampleIdx,:)),obj.data{simID-1},'uni',0);
            D.response = nan(size(D.response)); D = rmfield(D,{'feedbackType','RT','repeatNum'});
            D.stimulus = D.stimulus(:,1:2);
            
            switch(profile)
                case 1
                    %BO_NC50 baseline model every session has same
                    %parameters
                    %bias is affected in vis and m2 areas
                    %sens is also effected in vis
                    %no laser effects elsewhere
                    m = obj.getModel2('BO_NC50','all');
                    
                    bL = -2.2;
                    bR = -1.8;
                    sL = 1;
                    sR = 1;
                    n = 1.2;
                    c50 = 0.4;
                    
                    p_nL = [bL bR sL sR n n c50 c50];
                    p_L = {'Lvis', 1, [+1 -1  -0.5 0  0 0 0 0];...
                        'Rvis', 2, [-1 +1  0 -0.5  0 0 0 0];...
                        'Lm2',  5,  [+1 -1  0 0  0 0 0 0];...
                        'Rm2',  6,  [-1 +1  0 0  0 0 0 0]};
                    
                    %First simulate responses on trials without any laser
                    %effects
                    noLasIdx = ~any(D.areaIdx==[p_L{:,2}],2);
                    
                    p_hat = obj.calculatePhat(m.ZL,m.ZR,[],p_nL,D.stimulus(noLasIdx,:));
                    [~,resp]=max(mnrnd(1,p_hat),[],2);
                    D.response(noLasIdx) = resp;
                    
                    
                    %Now go through the trials for each laser-effect site
                    %and create trials
                    
                    for area = 1:size(p_L,1)
                        LasIdx = D.areaIdx==p_L{area,2};
                        
                        p_hat = obj.calculatePhat(m.deltaZL,m.deltaZR,p_nL,p_L{area,3},D.stimulus(LasIdx,:));
                        [~,resp]=max(mnrnd(1,p_hat),[],2);
                        D.response(LasIdx) = resp;
                    end
                    obj.data{simID} = D;
                    
                    tab = tabulate(D.response);
                    tab = tab(:,3)/100;
                    obj.guess_bpt(simID)=sum(tab.*log2(tab));
                    
            end
            
            
            
            
            
        end
        
        function obj = simulateMechanistic(obj,params)
            if ~any(strcmp(obj.names,'<<SIMULATION>>')==1)
                simID = length(obj.names)+1;
                
                obj.names{simID}  = '<<SIMULATION>>';
                obj.data{simID} = struct;
                obj.guess_bpt(simID) = nan;
            else
                simID = length(obj.names);
            end
            
            D = obj.data{simID-1};
            D.areaIdx(D.laserIdx==0)=0;
            D = getrow(D,D.areaIdx<=6);
            numTrials = length(D.response);
            
            %Go through each trial, simulate mechanistic model, then
            %replace the response with the model response
            for t = 1:numTrials
                p_hat = simulateMechanisticFcn(D.stimulus(t,1:2),params,D.areaIdx(t),obj.bilateral_flag);
                D.response(t,1) = find(mnrnd(1,p_hat));
            end
            
            obj.data{simID} = D;
        end
        
        function obj = removeSimulation(obj)
            if any(strcmp(obj.names,'<<SIMULATION>>')==1)
                simID = length(obj.names);
                obj.data(simID)=[];
                obj.guess_bpt(simID) = [];
                obj.names(simID) = [];
            end
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
        
        function obj = nonLaserBootstrap(obj) %OBSOLETE
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
        %
        function obj = fitNull(obj) %TO UPDATE
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
            numIter = 300;
            
            obj.fitData.deltaParamsNULL = cell(1,numSubjects);
            figure('name','Sampling from null distribution of laser parameters');
            for n = 1:numSubjects
                %Get all non-laser trials
                tab = tabulate(obj.data{n}.laserIdx);
                proportion = median(tab(2:end,3)/100);
                
                subplot(1,numSubjects,n); title(obj.names{n}); hold on;
                
                numLP = size(obj.fitData.params{n},2);
                if numLP > 2
                    scatter(obj.fitData.deltaParams{n}(:,1),obj.fitData.deltaParams{n}(:,3),30,cols,'filled');
                elseif numLP == 2
                    scatter(obj.fitData.deltaParams{n}(:,1),obj.fitData.deltaParams{n}(:,2),30,cols,'filled');
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
                    sessions = unique(obj.data{n}.sessionID(nL));
                    p_nL = nan(length(sessions),length(obj.model.pLabels));
                    for s = 1:length(sessions)
                        cont = obj.data{n}.stimulus(nL & obj.data{n}.sessionID==sessions(s),1:2);
                        resp = obj.data{n}.response(nL & obj.data{n}.sessionID==sessions(s));
                        
                        p_nL(s,:) = obj.util_optim(obj.model.ZL,obj.model.ZR,[],cont,resp,obj.model.LB,obj.model.UB);
                    end
                    p_nL(:,obj.model.LB==0 & obj.model.UB==0) = p_nL(:,find(obj.model.LB==0 & obj.model.UB==0) - 1);
                    
                    LB = -inf(1,length(obj.model.pLabels));
                    UB = inf(1,length(obj.model.pLabels));
                    LB(~obj.model.deltaP)=0;
                    UB(~obj.model.deltaP)=0;
                    
                    cont = obj.data{n}.stimulus(q_idx,1:2);
                    resp = obj.data{n}.response(q_idx);
                    sess = obj.data{n}.sessionID(q_idx);
                    
                    %offset on parameters
                    offset = p_nL(sess,:);
                    p = obj.util_optim(obj.model.deltaZL,obj.model.deltaZR,offset,cont,resp,LB,UB);
                    obj.fitData.deltaParamsNULL{n}(iter,:) = p;
                    
                    if numLP > 2
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
            obj.guess_bpt(1:numSubjects-1)=[];
            obj.names(1:numSubjects-1)=[];
        end
        
        function obj = clearNull(obj)
            obj.significance_flag=0;
            try
                obj.fitData = rmfield(obj.fitData,{'paramsNull'});
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
        
        function obj = crossval_LaserEffect(obj,n)
            %Fits complete nonlaser set of data and only does crossval on
            %the trials with laser inactivation
            
            perSession = 0;
            
            
            obj.cvData = [];
            
            areas = {%'V1', [1 2];
                %'M2', [5 6];
                %'S1', [3 4];
                'all' (1:9)};
            
            InactivationSites = [areas{:,2}];
            InactivationSites = [1:52];
            
            baseModel = 'BO_NC50';
            laserModels = {'bias','sens','^n','^c50','bias+sens','bias+^n','bias+^c50','bias+sens+^c50','bias+sens+^n','bias+^n+^c50','bias+sens+^n+^c50'};
            
            %             numLas = 2^length(laserModels) - 1;
            numLas = length(laserModels);
            
            MODELSTRINGS = cell(1,numLas + 1);
            MODELSTRINGS{end}='noLaserTerms';
            
            if isempty(obj.cvData)
                obj.cvData.L = cell(1,length(obj.names));
                obj.cvData.loglik = cell(1,length(obj.names));
                obj.cvData.areas = areas;
            end
            
            img=imread('D:\kirkcaldie_brain_BW_outline.png');
            
            %First fit the whole nonlaser dataset
            NL = getrow(obj.data{n},obj.data{n}.areaIdx==999);
            m = obj.getModel2(baseModel,'');
            
            sessions = unique(NL.sessionID);
            p_nL = nan(length(sessions),length(m.pLabels));
            disp(['Non-laser model loglik on entire dataset']);
            if perSession==1
                for s = 1:length(sessions)
                    cont = NL.stimulus(NL.sessionID==s,1:2);
                    resp = NL.response(NL.sessionID==s);
                    
                    p_nL(s,:) = obj.util_optim(m.ZL,m.ZR,[],cont,resp,m.LB,m.UB);
                    
                    %Print loglik
                    loglik=obj.calculateLogLik(p_nL(s,:),m.ZL,m.ZR,[],cont,resp)-obj.guess_bpt(n);
                    disp(['Session ' num2str(s) ': ' num2str(loglik)]);
                end
            else
                cont = NL.stimulus(:,1:2);
                resp = NL.response;
                p_nL = obj.util_optim(m.ZL,m.ZR,[],cont,resp,m.LB,m.UB);
                loglik=obj.calculateLogLik(p_nL,m.ZL,m.ZR,[],cont,resp)-obj.guess_bpt(n);
                disp(['Pooled-session: ' num2str(loglik)]);
                p_nL = repmat(p_nL,length(sessions),1);
            end
            p_nL(:,m.LB==0 & m.UB==0) = p_nL(:,find(m.LB==0 & m.UB==0) - 1);
            
            
            %Now do crossval on the laser part of the data
            numFolds = 10;
            %             L = getrow(obj.data{n},any(obj.data{n}.areaIdx == InactivationSites,2));
            L = getrow(obj.data{n},any(obj.data{n}.laserIdx == InactivationSites,2));
            L.p_hat = nan(size(L.response,1),numLas+1);
            L.p_hat_trainingset = nan(size(L.p_hat));
            
            %             cv = cvpartition(L.areaIdx,'KFold',numFolds);
            cv = cvpartition(L.laserIdx,'KFold',numFolds);
            
            loglik_perfold = [];
            for fold = 1:cv.NumTestSets
                disp(['fold: ' num2str(fold) '/' num2str(numFolds)]);
                for site=1:length(InactivationSites)
                    %                     trainIdx = cv.training(fold) & L.areaIdx==InactivationSites(site);
                    %                     testIdx = cv.test(fold) & L.areaIdx==InactivationSites(site);
                    
                    trainIdx = cv.training(fold) & L.laserIdx==InactivationSites(site);
                    testIdx = cv.test(fold) & L.laserIdx==InactivationSites(site);
                    
                    trainC = L.stimulus(trainIdx,1:2);
                    testC = L.stimulus(testIdx,1:2);
                    trainR = L.response(trainIdx);
                    testR = L.response(testIdx);
                    trainSID = L.sessionID(trainIdx);
                    testSID = L.sessionID(testIdx);
                    
                    %Fit on training laser data
                    trainOffset = p_nL(trainSID,:);
                    ph = obj.calculatePhat(m.deltaZL,m.deltaZR,trainOffset,zeros(1,length(m.pLabels)),trainC);
                    L.p_hat_trainingset(trainIdx,numLas + 1) = ph(:,1).*(trainR==1) + ph(:,2).*(trainR==2) + ph(:,3).*(trainR==3);
                    
                    %Test model without any laser terms
                    testOffset = p_nL(testSID,:);
                    ph = obj.calculatePhat(m.deltaZL,m.deltaZR,testOffset,zeros(1,length(m.pLabels)),testC);
                    L.p_hat(testIdx,numLas + 1) = ph(:,1).*(testR==1) + ph(:,2).*(testR==2) + ph(:,3).*(testR==3);
                    
                    for laserModel = 1:numLas
                        %                         lasIdx = find(decimalToBinaryVector(laserModel,length(laserModels)));
                        lasIdx = laserModel;
                        
                        laserString = strjoin(laserModels(lasIdx),'+');
                        MODELSTRINGS{1,laserModel} = [laserString];
                        m = obj.getModel2(baseModel,laserString);
                        LB = -inf(1,length(m.pLabels));
                        UB = inf(1,length(m.pLabels));
                        LB(~m.deltaP)=0;
                        UB(~m.deltaP)=0;
                        
                        %Fit laser model
                        p = obj.util_optim(m.deltaZL,m.deltaZR,trainOffset,trainC,trainR,LB,UB);
                        ph = obj.calculatePhat(m.deltaZL,m.deltaZR,trainOffset,p,trainC);
                        L.p_hat_trainingset(trainIdx,laserModel) = ph(:,1).*(trainR==1) + ph(:,2).*(trainR==2) + ph(:,3).*(trainR==3);
                        
                        %Test laser model
                        ph = obj.calculatePhat(m.deltaZL,m.deltaZR,testOffset,p,testC);
                        L.p_hat(testIdx,laserModel) = ph(:,1).*(testR==1) + ph(:,2).*(testR==2) + ph(:,3).*(testR==3);
                        %
                    end
                    
                end
                
                
                loglik = mean(log2(L.p_hat(cv.test(fold),:)),1);
                loglik_perfold(fold,:) = loglik(1:end-1) - loglik(end);
                
            end
            
            obj.cvData.L{n} = L;
            
            %             figure('name',obj.names{n},'color','w')
            %             for a = 1:size(areas,1)
            %
            %                 idx = any(L.areaIdx==areas{a,2},2);
            %
            %                 loglik = mean(log2(L.p_hat(idx,:)));
            %
            %                 obj.cvData.loglik{n}(a,:) = loglik;
            %
            %                 h(a)=subplot(1,size(areas,1),a);
            %                 ax = bar(loglik(1:end-1));
            %                 ax.BaseValue = loglik(end);
            %                 ax.EdgeAlpha=0; ax.BarWidth=0.95;
            %
            %                 set(gca,'XTickLabel',MODELSTRINGS(1:end-1),'XTickLabelRotation',45);
            %                 title(areas{a,1});
            %
            %                 if a == 1
            %                     ylabel({'log likelihood','baseline = value without laser terms'});
            %                 end
            %
            %                 set(gca,'box','off');
            %
            %             end
            %             linkaxes(h,'y');
            %
            obj.cvData.MODELSTRINGS = MODELSTRINGS;
            
            
            
            figure('name',obj.names{n},'color','w')
            %             for a = 1:size(areas,1)
            
            loglik = mean(loglik_perfold,1);
            stderr = std(loglik_perfold)/sqrt( size(loglik_perfold,1) );
            
            
            %                 loglik = mean(log2(L.p_hat),1);
            %                 laserloglik = loglik(1:end-1);
            
            [sorted_laserloglik,sorted_idx] = sort(loglik);
            
            %                 h(a)=subplot(1,size(areas,1),a);
            ax = barh(sorted_laserloglik);
            %                 ax.BaseValue = nonlaserloglik;
            ax.EdgeAlpha=0; ax.BarWidth=0.9;
            ax.FaceColor=[1 1 1]*0.9;
            
            
            set(gca,'YTickLabel',MODELSTRINGS(sorted_idx),'YTickLabelRotation',0);
            %                 title(areas{a,1});
            
            %                 if a == 1
            xlabel({'log likelihood','relative to model without perturbations'});
            %                 end
            
            set(gca,'box','off','ydir','reverse');
            
            
            hold on;
            er = errorbar(loglik(sorted_idx),1:length(loglik),stderr(sorted_idx),'horizontal');
            er.CapSize=0;
            er.LineStyle='none';
            er.LineWidth=3;
            er.Color=[1 1 1]*0.3;
            %             end
            %             linkaxes(h,'y');
            
            
        end
        
        function plotlasercrossval(obj,id1,id2)
            
            if isempty(obj.cvData)
                error('Do crossval on all data');
            end
            
            %             subjs = find(~cellfun(@isempty,obj.cvData.loglik));
            %             for n = 1:length(subjs)
            %                 figure('name',obj.names{subjs(n)},'color','w');
            %
            %                 loglik = obj.cvData.loglik{subjs(n)};
            %
            %                 for a = 1:size(obj.cvData.areas,1)
            %                     h(a)=subplot(1,size(obj.cvData.areas,1),a);
            %                     ax = bar(loglik(a,1:end-1));
            %                     ax.BaseValue = loglik(a,end);
            %                     ax.EdgeAlpha=0; ax.BarWidth=0.95;
            %
            %                     set(gca,'XTickLabel',obj.cvData.MODELSTRINGS(1:end-1),'XTickLabelRotation',45);
            %                     title(obj.cvData.areas{a,1});
            %
            %                     if a == 1
            %                         ylabel({'log likelihood','baseline = value without laser terms'});
            %                     end
            %
            %                 end
            %                 linkaxes(h,'y');
            %             end
            %
            %
            %
            
            deltaLL = [];
            
            
            
            for n = 1:(length(obj.names)-1)
                deltaLL(n,:) = (obj.cvData.loglik{n}(:,id1) - obj.cvData.loglik{n}(:,id2))';
            end
            
            figure('color','w');
            
            h = scatterhist(deltaLL(:,1),deltaLL(:,2),'Group',obj.names(1:end-1)','Marker','.','MarkerSize',50);
            set(h(1),'box','off');
            hold on;
            
            boxplot(h(2),deltaLL(:,1),'orientation','horizontal');
            boxplot(h(3),deltaLL(:,2),'orientation','horizontal'); view(h(3),[270,90]);
            set(h(2:3),'XTickLabel','','YTickLabel','','box','on');
            
            fimplicit(h(1),@(x,y)(x-y))
            %             errorbar(mean(deltaLL),1.96*std(deltaLL)/sqrt(size(deltaLL,1)));
            xlabel('V1');
            ylabel('M2');
            
            
            title({[obj.cvData.MODELSTRINGS{id1} ' relative to ' obj.cvData.MODELSTRINGS{id2}],...
                ['H0: V1=0 : p=' num2str( signrank(deltaLL(:,1)) )],...
                ['H0: M2=0 : p=' num2str( signrank(deltaLL(:,2)) )],...
                ['H0: V1=M2 : p=' num2str( signrank(deltaLL(:,1),deltaLL(:,2)) )]});
        end
        
        
        function plotFit_NonLaserOneSession(obj,n,sess)
            if isempty(sess)
                p_nL = median(obj.fitData.params{n},1);
                idx = obj.data{n}.laserIdx==0;
                if obj.fitData.perSession==1
                    warning('Model fits cannot be visualised well here as this requires averaging nonlaser params');
                end
            else
                p_nL = obj.fitData.params{n}(sess,:);
                idx = obj.data{n}.laserIdx==0 & obj.data{n}.sessionID==sess;
            end
            
            %             [ZL_nL,ZR_nL,LB,~] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
            
            cont = obj.data{n}.stimulus(idx,1:2);
            resp = obj.data{n}.response(idx);
            
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
                    p_pred(cl,cr,:) = obj.calculatePhat(obj.model.ZL,obj.model.ZR,[],p_nL,[cVals(cl) cVals(cr)]);
                end
            end
            
            %             labels = {'pL','pR','pNG'};
            %             for r = 1:3
            %                 subplot(3,4,r);
            %                 imagesc(p(:,:,r)); set(gca,'ydir','normal'); title(['Actual ' labels{r}]);
            %                 caxis([0 1]);
            %                 set(gca,'xtick',1:length(cVals),'ytick',1:length(cVals),'xticklabels',cVals,'yticklabels',cVals);
            %                 xlabel('CR'); ylabel('CL');
            %                 subplot(3,4,r+4);
            %                 imagesc(p_pred(:,:,r)); set(gca,'ydir','normal'); title(['Pred ' labels{r}]);
            %                 caxis([0 1]);
            %                 set(gca,'xtick',1:length(cVals),'ytick',1:length(cVals),'xticklabels',cVals,'yticklabels',cVals);
            %             end
            %             axis(get(gcf,'children'),'square');
            %             set(get(gcf,'children'),'ytick',cVals,'xtick',cVals);
            
            
            %pych curve: pedestal representation
            
            
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
                        %                         l(2)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,1));
                        %                         l(3)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,2));
                        set(l,'Color',cols(ch,:),'Linewidth',0.5);
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
                
                
                title(['pedestal= ' num2str(cVals(ped))]);
                axis square;
                if ped == numPedestals
                    ylabel('Choice probability');
                    xlabel('\Delta contrast from pedestal');
                end
                if ped < numPedestals
                    set(gca,'xtick','','xcolor','w');
                end
                
            end
            
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
        
        function plot_history(obj,type,n)
            d = obj.data{n};
            d = getrow(d,d.laserType==0); %Only non-laser trials
            
            cols = [        0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

           
            
            choiceLabels = {'Left','Right','NoGo'};
            
            switch(type)
                case 'stim'
                    e = {};
                    e{1} = getrow(d,d.prev_stimulus(:,1) > d.prev_stimulus(:,2));
                    e{2} = getrow(d,d.prev_stimulus(:,1) < d.prev_stimulus(:,2));
                    e{3} = getrow(d,d.prev_stimulus(:,1) == 0 & d.prev_stimulus(:,2) == 0);
                    
                    legendLabel = {'Stimulus CL>CR','Stimulus CR>CL','Stimulus Zero'};
                case 'choice'
                    e = {};
                    e{1} = getrow(d,d.prev_response == 1);
                    e{2} = getrow(d,d.prev_response == 2);
                    e{3} = getrow(d,d.prev_response == 3);
                    
                    legendLabel = {'Response Left','Response Right','Response NoGo'};
                case 'reward'
                    e = {};
                    e{1} = getrow(d,d.prev_feedback == 1);
                    e{2} = getrow(d,d.prev_feedback == -1);
                    
                    legendLabel = {'Rewarded','Punished'};
                case 'rewardGo'
                    e = {};
                    e{1} = getrow(d,d.prev_response < 3 & d.prev_feedback == 1);
                    e{2} = getrow(d,d.prev_response < 3 & d.prev_feedback == -1);
                    
                    legendLabel = {'Rewarded Go','Punished Go'};
                case 'rewardLeft'
                    e = {};
                    e{1} = getrow(d,d.prev_response == 1 & d.prev_feedback == 1);
                    e{2} = getrow(d,d.prev_response == 1 & d.prev_feedback == -1);
                    
                    legendLabel = {'Rewarded Left','Punished Left'};
                    
                case 'rewardRight'
                    e = {};
                    e{1} = getrow(d,d.prev_response == 2 & d.prev_feedback == 1);
                    e{2} = getrow(d,d.prev_response == 2 & d.prev_feedback == -1);
                    
                    legendLabel = {'Rewarded Right','Punished Right'};
            end
            
            g = {};
            for subset = 1:length(e)
                g{subset} = GLM(e{subset}).setModel('C50-subset').fit;
            end
            
            %Overlay fitted curves
            cVal = unique(d.stimulus(:));
            numPedestals = length(cVal) - 1;
            
            %Plot non-laser predictions and data
            figure('color','w','name',obj.names{n});
            for ped = 1:numPedestals
                
                testCont = [linspace( 0.54 - cVal(ped) + 0.1 ,0,100)' zeros(100,1); zeros(100,1) linspace(0,0.54 - cVal(ped) + 0.1,100)'] + cVal(ped);
                
                phat = {};
                for subset = 1:length(e)
                    phat{subset} = g{subset}.calculatePhat(g{subset}.parameterFits,testCont);
                end
                
                ii=1;
                for r = [1 3 2]
                    subplot(numPedestals,3,3*ped -3 + ii); hold on;
                    
                    for subset = 1:length(e)
                        hx=plot(diff(testCont,[],2),phat{subset}(:,r),'linewidth',1.5);
                        hx.Color = cols(subset,:);
                    end

                    xlim([-1 1]*(max(cVal)+0.1)); ylim([0 1]);
                    
                    ii=ii+1;
                    
                    if ped==1
                        title(choiceLabels{r});
                    end
                    
                    if ped < numPedestals
                        set(gca,'xtick','','xcolor','w');
                    end
                    
                    if r > 1
                        set(gca,'ytick','','ycolor','w');
                    end
                end
                
                
            end
            
            legend(legendLabel);
            legend boxoff;
            
            
            
      
            
        end
        
        function plotFit_NonLaser(obj)
            numSubjects = length(obj.names);
            
            %Plot session-by-session bias and sensitivity parameters
            figure('name','Session by session non-laser bias and sensitivity changes','color','w');
            i=1;
            for p = [1 2 3 4 5 7]
                h(i)=subplot(3,2,i); hold on;
                
                for n = 1:(length(obj.names)-1)
                    numSessions = length(obj.expRefs{n});
                    %                     ax=notBoxPlot(obj.fitData.params{n}(:,p),n,'jitter',0.5,'markMedian',true);
                    %                     ax.data.MarkerSize=3;
                    %                     ax.data.MarkerFaceColor=[1 1 1]*0.2;
                    %                     ax.data.MarkerEdgeColor=[1 1 1]*0.2;
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
            
            
            %Plot non-laser predictions against actual psychometric data.
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
                    idx = obj.data{n}.laserIdx==0 & obj.data{n}.sessionID==s;
                    
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
                                %                                 l(2)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,1));
                                %                                 l(3)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,2));
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
        
        
        function plotFit_Laser(obj,varargin)
            if isempty(varargin)
                subjID = 1:length(obj.names);
            else
                subjID = varargin{1};
            end
            numSubjects = length(obj.names);
            
            img=imread('D:\kirkcaldie_brain_BW_outline.png');
            cols = {[0 0.4470 0.7410],...
                [0.8500 0.3250 0.0980],...
                [0.4940    0.1840    0.5560]};
            
            if obj.fitData.groupAreas==0 && length(subjID)>1
                %Plot laser parameter maps
                figure('name','Laser parameter maps','color','w');
                for n = 1:length(obj.names)
                    p_L = obj.fitData.deltaParams{n};
                    
                    P_old_mag = abs(mean(obj.fitData.params{n},1));
                    plotVal = cell(1,size(p_L,2));
                    for p = 1:size(p_L,2)
                        plotVal{p} = p_L(:,p);
                        % %                     plotVal{p} = plotVal{p}/std(obj.fitData.paramsNull{n}(:,p));
                        %                     plotVal{p} = (P_old_mag(p) + plotVal{p})./P_old_mag(p);
                        
                    end
                    %                 keyboard;
                    
                    param_labels = obj.model.pLabels;
                    
                    if obj.significance_flag == 1
                        pvalue = obj.sig{n};
                        dotSize = 35 + (pvalue<0.05)*80 + (pvalue<0.01)*100;
                    else
                        dotSize = ones(size(p_L))*150;
                    end
                    
                    for i = 1:size(plotVal,2)
                        subplot(numSubjects,size(plotVal,2),size(plotVal,2)*n-size(plotVal,2) + i);
                        
                        ap = obj.inactivationCoords(:,2);
                        ml = obj.inactivationCoords(:,1);
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
                            title(['\Delta ' param_labels{i}]);
                        end
                    end
                    
                    cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                        linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                    colormap(flipud(cmap));
                end
            end
            
            %             figure('color','w','name','ALLSUBJ');
            %             for p = 1:size(plotVal,2)
            %                 subplot(size(plotVal,2)/2,2,p);
            %                 imX=imagesc(linspace(-4.5,4.5,100),linspace(3.75,-5.2,100),img);
            %                 set(gca,'ydir','normal');
            %                         set(imX,'alphadata',0.7); hold on;
            %
            %                 scatter(ml,ap,dot,plotVal{p},'s','filled'); axis equal;
            %                 caxis([-1 1]*max(abs(plotVal{p}))); colorbar;
            %                 set(gca,'box','off','ytick','','xtick','','xcolor','w','ycolor','w');
            %                 ylim([-6,4]);
            %             end
            %             cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
            %                         linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
            %                     colormap(flipud(cmap));
            
            if obj.fitData.groupAreas==0 && length(subjID)==1
                figure('name',['Laser parameter maps' obj.names{subjID}],'color','w');
                p_L = obj.fitData.deltaParams{subjID};
                param_labels = obj.model.pLabels;
                plotVal = cell(1,size(p_L,2));
                for p = 1:size(p_L,2)
                    plotVal{p} = p_L(:,p);
                end
                dotSize = ones(size(p_L))*150;
                for i = 1:size(plotVal,2)
                    subplot(size(plotVal,2)/2,2,i);
                    
                    ap = obj.inactivationCoords(:,2);
                    ml = obj.inactivationCoords(:,1);
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
                    
                    title(['\Delta ' param_labels{i}]);
                end
                
                cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
                colormap(flipud(cmap));
            end
            
            figure('color','w');
            params = find(obj.model.deltaP);
            dP = cat(3,obj.fitData.deltaParams{:});
            for p = 1:length(params)
                subplot(length(params),1,p);
                title(obj.model.pLabels{params(p)});
                notBoxPlot(squeeze(dP(:,params(p),1:numSubjects-1))');
                
                hold on;
                line(get(gca,'XLim'),[0 0]);
                if p < length(params)
                    set(gca,'xtick','','xcolor','w');
                end
            end
            if obj.fitData.groupAreas==1
                set(gca,'XTickLabel',{'Left V1','Right V1','Left S1','Right S1','Left M2','Right M2'});
            end
            
            %             %Plot psych curves showing the laser effect at a few sites, for
            %             %each mouse separately
            if obj.fitData.perSession==1 && obj.fitData.groupAreas == 0
                error('psych-curve method not good for per-session fits because I have to average nonlaser params over sessions');
            end
            %
            if obj.bilateral_flag==0
                siteLabels = {'Right V1','Right M2','Right Barrel'};
                if obj.fitData.groupAreas==1
                    siteID = [2 6 4];
                else
                    siteID = [10 15 48 49 ];%25 32];
                end
            elseif obj.bilateral_flag==1
                siteLabels = {'V1','M2','Barrel'};
                if obj.fitData.groupAreas==1
                    siteID = [2 6 4];
                else
                    siteID = [7 24 ];
                end
            end
            
            choiceLabels = {'Choose Left','Choose Right','NoGo'};
            
            for n = subjID
                %                 figure('name',[obj.names{n} ' ' obj.model.name{1} ' [' obj.model.name{2} ']'],'color','w');
                
                
                for site = 1:length(siteLabels)
                    
                    %                     subplot(1,length(siteLabels),site);
                    
                    cont = obj.data{n}.stimulus(:,1:2);
                    resp = obj.data{n}.response;
                    cVals = unique(unique(cont(:,1:2)));
                    numPedestals = length(cVals);
                    
                    p_nL = mean(obj.fitData.params{n});
                    nL_idx = obj.data{n}.laserIdx==0;
                    
                    if obj.fitData.groupAreas==1
                        p_L = obj.fitData.deltaParams{n}(siteID(site),:);
                        L_idx = obj.data{n}.areaIdx==siteID(site);
                    else
                        p_L = obj.fitData.deltaParams{n}(siteID(site),:);
                        L_idx = obj.data{n}.laserIdx==siteID(site);
                    end
                    
                    figure('name',[obj.names{n} ' ' siteLabels{site} ' ' obj.model.name{1} ' [' obj.model.name{2} ']'],'color','w');
                    %Plot non-laser predictions and data
                    for ped = 1:numPedestals
                        cnt = cont(nL_idx,:);
                        res = resp(nL_idx);
                        
                        uC = unique(diff(cnt(min(cnt,[],2)==cVals(ped),:),[],2));
                        testCont = [linspace(max(abs(uC))+0.1,0,100)' zeros(100,1); zeros(100,1) linspace(0,max(abs(uC))+0.1,100)'] + cVals(ped);
                        p_hat = obj.calculatePhat(obj.model.ZL,obj.model.ZR,[],p_nL,testCont);
                        %                         set(gca,'ColorOrderIndex',1);
                        %                         plot(diff(testCont,[],2),p_hat,'linewidth',0.5);
                        
                        ii=1;
                        for r = [1 3 2]
                            subplot(numPedestals,3,3*ped -3 + ii); hold on;
                            plot(diff(testCont,[],2),p_hat(:,r),'k','linewidth',0.5);
                            xlim([-1 1]*(max(cVals)+0.1)); ylim([0 1]);
                            
                            ii=ii+1;
                            
                            if ped==1
                                title(choiceLabels{r});
                            end
                        end
                        
                        %Plot pooled data as well
                        ped_idx = min(cnt,[],2)==cVals(ped);
                        ped_c_diff = diff(cnt(ped_idx,:),[],2);
                        ped_r = res(ped_idx);
                        uC = unique(ped_c_diff);
                        ph=[];
                        for c = 1:length(uC)
                            r = ped_r(ped_c_diff==uC(c));
                            [ph(c,:),pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
                            
                            ii = 1;
                            for r=[1 3 2]
                                subplot(numPedestals,3,3*ped -3 + ii); hold on;
                                l=line([1 1]*uC(c),pci(r,:));
                                plot(uC(c),ph(c,r),'k.','markersize',10);
                                %                                 l(2)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(r,1));
                                %                                 l(3)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(r,2));
                                set(l,'Color','k','Linewidth',0.5);
                                %                             l.Color = cols{ch};
                                %                             l.LineWidth=1;
                                ii = ii+1;
                            end
                        end
                        
                        %
                        uC = unique(diff(cnt(min(cnt,[],2)==cVals(ped),:),[],2));
                        testCont = [linspace(max(abs(uC))+0.1,0,100)' zeros(100,1); zeros(100,1) linspace(0,max(abs(uC))+0.1,100)'] + cVals(ped);
                        %Overlay laser predictions
                        p_hat = obj.calculatePhat(obj.model.deltaZL,obj.model.deltaZR,p_nL,p_L,testCont);
                        %                         set(gca,'ColorOrderIndex',1);
                        %                         plot(diff(testCont,[],2),p_hat,'linewidth',3);
                        ii=1;
                        for r = [1 3 2]
                            subplot(numPedestals,3,3*ped -3 + ii); hold on;
                            h=plot(diff(testCont,[],2),p_hat(:,r),'r','linewidth',0.5);
                            %                             h.Color=siteCols{site};
                            ii=ii+1;
                        end
                        
                        
                        %And actual laser data
                        %                         keyboard
                        cnt = cont(L_idx,:);
                        res = resp(L_idx);
                        ped_idx = min(cnt,[],2)==cVals(ped);
                        ped_c_diff = diff(cnt(ped_idx,:),[],2);
                        ped_r = res(ped_idx);
                        uC = unique(ped_c_diff);
                        ph=[];
                        for c = 1:length(uC)
                            r = ped_r(ped_c_diff==uC(c));
                            [ph(c,:),pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
                            
                            ii=1;
                            for r=[1 3 2]
                                subplot(numPedestals,3,3*ped -3 + ii); hold on;
                                l=line([1 1]*uC(c),pci(r,:));
                                h=plot(uC(c),ph(c,r),'r.','markersize',10);
                                %                                 h.Color=siteCols{site};
                                %                                 l(2)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,1));
                                %                                 l(3)=line([uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,2));
                                set(l,'Color','r','Linewidth',0.5);
                                %                             l.Color = cols{ch};
                                %                             l.LineWidth=1;
                                ii = ii+1;
                                
                                if r>1
                                    set(gca,'ytick','','ycolor','w');
                                end
                                
                                if ped<numPedestals
                                    set(gca,'xtick','','xcolor','w');
                                end
                            end
                        end
                        
                    end
                    
                    %                     keyboard;
                end
                
            end
            
            %Plot change to C50 function, if the model specifies n
            %or c50 to change
            
            
            %             if any(sum(abs(obj.fitData.deltaParams{end}(:,5:end)),1)>0)
            % %                 keyboard;
            %                 figure('name',['CFN FCN CHANGE' obj.model.name{1} ],'color','w');
            %
            %                 ii=1;
            %                 for n = subjID
            %                     for site = 1:length(siteLabels)
            %                         subplot(length(subjID),length(siteLabels),length(siteLabels)*(ii-1) + site);
            %                         hold on;
            %
            %                         p_nL=obj.fitData.params{n}(1,5:end);
            %                         p_L=(p_nL).*2.^obj.fitData.deltaParams{n}(site,5:end);
            %
            %                         if length(p_nL)==4
            %                             h=ezplot(@(c)(cfn(c,p_nL(1),p_nL(3))),[0 1]); h.Color=[0 0 0];
            %                             h=ezplot(@(c)(cfn(c,p_L(1),p_L(3))),[0 1]); h.Color=[1 0 0];
            %
            %                             h=ezplot(@(c)(cfn(c,p_nL(2),p_nL(4))),[0 1]); h.Color=[0 0 0]; h.LineStyle='--';
            %                             h=ezplot(@(c)(cfn(c,p_L(2),p_L(4))),[0 1]); h.Color=[1 0 0];h.LineStyle='--';
            %                         else
            %                             h=ezplot(@(c)(obj.model.cfn(c,p_nL(1))),[0 1]); h.Color=[0 0 0];
            %                             h=ezplot(@(c)(obj.model.cfn(c,p_L(1))),[0 1]); h.Color=[1 0 0];
            %
            %                             h=ezplot(@(c)(obj.model.cfn(c,p_nL(2))),[0 1]); h.Color=[0 0 0];h.LineStyle='--';
            %                             h=ezplot(@(c)(obj.model.cfn(c,p_L(2))),[0 1]); h.Color=[1 0 0];h.LineStyle='--';
            %
            %                         end
            %
            %                         title(siteLabels{site}); axis square;
            %
            %                         if site==1
            %                             ylabel(obj.names{n});
            %                         end
            %                     end
            %
            %                     ii=ii+1;
            %                 end
            %
            %                 legend('f(cL)','f(cL) + laser','f(cR)','f(cR) + laser');
            %             end
            %
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
                %                 maxC = max(obj.data{n}.stimulus);
                %                 c_0 = sum(obj.data{n}.stimulus,2)==0;
                %                 c_L = obj.data{n}.stimulus(:,1) == maxC(1) & obj.data{n}.stimulus(:,2) == 0;
                %                 c_R = obj.data{n}.stimulus(:,2) == maxC(2) & obj.data{n}.stimulus(:,1) == 0;
                %                 c_LR = obj.data{n}.stimulus(:,1) == maxC(1) & obj.data{n}.stimulus(:,2) == maxC(2);
                %                 stim = c_0 + 2*c_L + 3*c_R + 4*c_LR;
                %                 stim_labels = {'C=[0 0]','High CL','High CR','High CL&CR'};
                
                % % %                 %quadrant conditions
                %                 c_0 = sum(obj.data{n}.stimulus,2)==0;
                %                 c_L = obj.data{n}.stimulus(:,1) > obj.data{n}.stimulus(:,2);
                %                 c_R = obj.data{n}.stimulus(:,2) > obj.data{n}.stimulus(:,1);
                %                 c_LR = ( obj.data{n}.stimulus(:,1) == obj.data{n}.stimulus(:,2) ) & ~c_0;
                %                 stim = c_0 + 2*c_L + 3*c_R + 4*c_LR;
                %                 stim_labels = {'C=[0 0]','CL > CR','CL < CR','CL = CR'};
                % %
                %Detection + CL=CR conditions
                c_0 = sum(obj.data{n}.stimulus,2)==0;
                c_L = obj.data{n}.stimulus(:,1) > 0 & obj.data{n}.stimulus(:,2)==0;
                c_R = obj.data{n}.stimulus(:,2) > 0 & obj.data{n}.stimulus(:,1)==0;
                c_LR = ( obj.data{n}.stimulus(:,1) == obj.data{n}.stimulus(:,2) ) & ~c_0;
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
                        cont = obj.data{n}.stimulus;
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
                
                ap = obj.inactivationCoords(:,2);
                ml = obj.inactivationCoords(:,1);
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
            %             ap = obj.inactivationCoords(:,2);
            %             ml = obj.inactivationCoords(:,1);
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
                cVal = unique(unique(D.stimulus(:,1:2)));
                p=nan(length(cVal),length(cVal),3);
                for cl = 1:length(cVal)
                    for cr = 1:length(cVal)
                        r = D.response(D.stimulus(:,1)==cVal(cl) & D.stimulus(:,2)==cVal(cr));
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
                [~,~,cValIdx]=unique(D.stimulus(:,1:2));
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
                areaList = {'left vis', 1; 'right vis', 2; 'left m2', 5; 'right m2', 6};
            else
                areaList = {'visual',2; 'm2',6};
            end
            
            %             if obj.nestedModel_flag == 1
            %                 xlab = 'log(pL/pR) given go';
            %                 ylab = 'log(pGO/pNG)';
            %             else
            xlab = 'log(pR/pNG)';
            ylab = 'log(pL/pNG)';
            %             end
            
            
            figure('color','w');
            for area=1:size(areaList,1)
                laserIdxs = areaList{area,2};
                cVal = unique(unique(D.stimulus(:,1:2)));
                
                p=nan(length(cVal),length(cVal),3);
                pLas=nan(length(cVal),length(cVal),3);
                C=nan(length(cVal),length(cVal),2);
                ZL_ci=nan(length(cVal),length(cVal),2);
                ZR_ci=nan(length(cVal),length(cVal),2);
                for cl = 1:length(cVal)
                    for cr = 1:length(cVal)
                        C(cl,cr,:) = cVal([cl cr]);
                        r = D.response(D.stimulus(:,1)==cVal(cl) & D.stimulus(:,2)==cVal(cr) & D.laserIdx==0);
                        p(cl,cr,:) = sum([r==1 r==2 r==3],1)/length(r);
                        
                        psim = squeeze(p(cl,cr,:));
                        n = length(r);
                        psimDraws=mnrnd(n,psim,5000)/n;
                        
                        ZL_ci(cl,cr,:) = quantile(log(psimDraws(:,1)./psimDraws(:,3)),[0.025 0.975]);
                        ZR_ci(cl,cr,:) = quantile(log(psimDraws(:,2)./psimDraws(:,3)),[0.025 0.975]);
                        
                        ZGO_ci(cl,cr,:) = quantile(log((1-psimDraws(:,3))./psimDraws(:,3)),[0.025 0.975]);
                        psimDraws_LR = bsxfun(@rdivide,psimDraws(:,1:2),sum(psimDraws(:,1:2),2));
                        ZLRG_ci(cl,cr,:) = quantile(log(psimDraws_LR(:,1)./psimDraws_LR(:,2)),[0.025 0.975]);
                        
                        r = D.response(D.stimulus(:,1)==cVal(cl) & D.stimulus(:,2)==cVal(cr) & arrayfun(@(l)(any(l==laserIdxs)),D.areaIdx));
                        %                         r = D.response(D.stimulus(:,1)==cVal(cl) & D.stimulus(:,2)==cVal(cr) & arrayfun(@(l)(any(l==laserIdxs)),D.laserIdx));
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
                
                %plot in log L/NG and log R/NG space, irrespective of
                %whether the model was fit in that space or the nested
                %logit space
                
                %                 if obj.nestedModel_flag == 0
                ZL = log(p(:,:,1)./p(:,:,3));
                ZR = log(p(:,:,2)./p(:,:,3));
                
                ZL_ciLOW = ZL_ci(:,:,1);
                ZL_ciHIGH = ZL_ci(:,:,2);
                ZR_ciLOW = ZR_ci(:,:,1);
                ZR_ciHIGH = ZR_ci(:,:,2);
                
                ZLLas = log(pLas(:,:,1)./pLas(:,:,3)); ZLLas(isinf(ZLLas))=NaN;
                ZRLas = log(pLas(:,:,2)./pLas(:,:,3)); ZRLas(isinf(ZRLas))=NaN;
                
                scatter(ZR(:),ZL(:),50,'o','filled'); hold on;
                plot(ZR,ZL,'k:');
                plot(ZR',ZL','k:');
                quiver(ZR(:),ZL(:),ZRLas(:)-ZR(:),ZLLas(:)-ZL(:),0);
                try
                    ezplot('2*exp(x)=1+exp(x)+exp(y)');
                    ezplot('2*exp(y)=1+exp(x)+exp(y)');
                    ezplot('2=1+exp(x)+exp(y)');
                catch
                    warning('error on ezplot of contours');
                end
                
                %                 pL = p(:,:,1);
                %                 pR = p(:,:,2);
                %                 pLLas = pLas(:,:,1);
                %                 pRLas = pLas(:,:,2);
                %                 scatter(pR(:),pL(:),50,'o','filled'); hold on;
                %                 quiver(pR(:),pL(:),pRLas(:)-pR(:),pLLas(:)-pL(:),0);
                %                 plot(pR,pL,'k:');
                %                 plot(pR',pL','k:');
                
                
                %                     l=line([ZR(:)'; ZR(:)'],[ZL_ciLOW(:)';ZL_ciHIGH(:)']);
                %                     set(l,'color',[0 0.4470 0.7410],'linestyle','-');
                %                     l=line([ZR_ciLOW(:)';ZR_ciHIGH(:)'],[ZL(:)'; ZL(:)']);
                %                     set(l,'color',[0 0.4470 0.7410],'linestyle','-');
                %                 else
                %                     pG(:,:,1) = p(:,:,1)./sum(p(:,:,1:2),3);
                %                     pG(:,:,2) = p(:,:,2)./sum(p(:,:,1:2),3);
                %                     pGLas(:,:,1) = pLas(:,:,1)./sum(pLas(:,:,1:2),3);
                %                     pGLas(:,:,2) = pLas(:,:,2)./sum(pLas(:,:,1:2),3);
                %
                %                     ZGO = log((1-p(:,:,3))./p(:,:,3));
                %                     ZLRG = log(pG(:,:,1)./pG(:,:,2));
                %
                %                     ZGO_ciLOW = ZGO_ci(:,:,1);
                %                     ZGO_ciHIGH = ZGO_ci(:,:,2);
                %                     ZLRG_ciLOW = ZLRG_ci(:,:,1);
                %                     ZLRG_ciHIGH = ZLRG_ci(:,:,2);
                %
                %                     ZGOLas = log((1-pLas(:,:,3))./pLas(:,:,3)); ZGOLas(isinf(ZGOLas))=NaN;
                %                     ZRLGLas = log(pGLas(:,:,1)./pGLas(:,:,2)); ZRLGLas(isinf(ZRLGLas))=NaN;
                %                     scatter(ZLRG(:),ZGO(:),50,'o','filled'); hold on;
                %                     plot(ZLRG([1:4 8 7 6 5 9:12 16 15 14 13]),ZGO([1:4 8 7 6 5 9:12 16 15 14 13]),'k:');
                %                     plot(ZLRG([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),ZGO([1 5 9 13 14 10 6 2 3 7 11 15 16 12 8 4]),'k:');
                %                     quiver(ZLRG(:),ZGO(:),ZRLGLas(:)-ZLRG(:),ZGOLas(:)-ZGO(:),0);
                %
                %                     l=line([ZLRG(:)'; ZLRG(:)'],[ZGO_ciLOW(:)';ZGO_ciHIGH(:)']);
                %                     set(l,'color',[0 0.4470 0.7410],'linestyle','-');
                %                     l=line([ZLRG_ciLOW(:)';ZLRG_ciHIGH(:)'],[ZGO(:)'; ZGO(:)']);
                %                     set(l,'color',[0 0.4470 0.7410],'linestyle','-');
                %
                %                     axis equal;
                %                 end
                axis equal;
                xlabel(['empirical ' xlab]); ylabel(['empirical ' ylab]); title(areaList{area,1});
                
                %overlay black dot for random nonlaser behaviour
                r = D.response(D.laserIdx==0);
                r=sum([r==1 r==2 r==3],1)/length(r);
                hold on;
                h = plot(log(r(2)/r(3)),log(r(1)/r(3)),'k.');
                hold off;
                h.MarkerSize=30;
                
                
                %                 try
                %                 if cell2mat(strfind(varargin,'withField')) == 1
                %                     cSet = {cVal,linspace(0,0.54,30)'};
                %                 else
                cSet = {cVal};
                %                 end
                
                h2(area)=subplot(2,size(areaList,1),area+size(areaList,1)); hold on;
                %                     if obj.nestedModel_flag == 0
                try
                    ezplot('2*exp(x)=1+exp(x)+exp(y)');
                    ezplot('2*exp(y)=1+exp(x)+exp(y)');
                    ezplot('2=1+exp(x)+exp(y)');
                catch
                    warning('error on ezplot of contours');
                end
                %                     end
                pNL = median(obj.fitData.params{subj},1);
                pL = median(obj.fitData.deltaParams{subj}(unique(obj.data{subj}.laserIdx(obj.data{subj}.areaIdx==laserIdxs)),:),1);
                
                for i = 1:length(cSet)
                    [cr,cl]=meshgrid(cSet{i});
                    cr = cr(:); cl = cl(:);
                    
                    
                    p_nonLaser = obj.calculatePhat(obj.model.ZL,obj.model.ZR,[],pNL,[cl cr]);
                    p_Laser = obj.calculatePhat(obj.model.deltaZL,obj.model.deltaZR,pNL,pL,[cl cr]);
                    
                    ZL = log(p_nonLaser(:,1)./p_nonLaser(:,3)); ZL = reshape(ZL,length(cVal),length(cVal));
                    ZR = log(p_nonLaser(:,2)./p_nonLaser(:,3)); ZR = reshape(ZR,length(cVal),length(cVal));
                    
                    ZLLas = log(p_Laser(:,1)./p_Laser(:,3)); ZLLas = reshape(ZLLas,length(cVal),length(cVal));
                    ZRLas = log(p_Laser(:,2)./p_Laser(:,3)); ZRLas = reshape(ZRLas,length(cVal),length(cVal));
                    %
                    %
                    %                         [ZLfcn,ZRfcn,LB] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,subj);
                    %                         ZL = ZLfcn([],pNL,[cl cr]); %ZL = reshape(ZL,length(cVal),length(cVal));
                    %                         ZR = ZRfcn([],pNL,[cl cr]); %ZR = reshape(ZR,length(cVal),length(cVal));
                    %
                    %                         [ZLfcn,ZRfcn] = obj.getModel(obj.fitData.laserModel,obj.fitData.biasMode,subj);
                    %
                    % %                         offset = repmat(pNL,length(cl),1);
                    %                         ZLLas = ZLfcn(pNL,pL,[cl cr]); %ZLLas = reshape(ZLLas,length(cVal),length(cVal));
                    %                         ZRLas = ZRfcn(pNL,pL,[cl cr]); %ZRLas = reshape(ZRLas,length(cVal),length(cVal));
                    %
                    if i ==1
                        scatter(ZR(:),ZL(:),50,'o','filled');
                        plot(ZR,ZL,'k:');
                        plot(ZR',ZL','k:');
                        quiver(ZR(:),ZL(:),ZRLas(:)-ZR(:),ZLLas(:)-ZL(:),0);
                    else
                        quiver(ZR(:),ZL(:),ZRLas(:)-ZR(:),ZLLas(:)-ZL(:));
                    end
                    
                end
                axis equal;
                xlabel(['predicted ' xlab]); ylabel(['predicted ' ylab]); title('');
                %                 catch
                %                 end
                
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
        
        function plotEffectOfTime(obj,n,type)
            d = obj.data{n};
            d = getrow(d, (d.stimulus(:,1)>0 | d.stimulus(:,2)>0) & d.stimulus(:,1)~=d.stimulus(:,2) );
            
            %Plots the effect of inactivating at different timepoints
            if length(unique(d.laserPower)) > 2
                powers = unique(d.laserPower);
                powers = powers(2:end);
                
                fprintf('<strong> Laser powers available, please select one </strong>\n');
                for i = 1:length(powers)
                    fprintf('%i) %0.1f mW\n',i,powers(i));
                end
                
                p = input('ID=');
                
                if p <= length(powers)
                    sessions = d.sessionID(d.laserPower==powers(p));
                    sessions = unique(sessions);
                    d = getrow(d, any(d.sessionID==sessions',2)); %Only get data from sessions using that power
                else
                    error('not valid power option');
                end
            end
            
            switch(type)
                case 'relStim'
                    d.plotDelay = d.laserOnset;
                    xlab = 'time of laser onset relative to stimulus onset';
                case 'relMove'
                    d = getrow(d,d.response<3);
                    d.plotDelay = d.RT - d.laserOnset;
                    xlab = 'time of movement onset relative to laser onset';
            end
            
            %Create figs
            fig = figure('name',obj.names{n},'color','w');
            ax_delays = subplot(5,1,1); hold on; title(ax_delays,'delays')
            ax_rtCorr = subplot(5,1,2); hold on; title(ax_rtCorr,'RT on correct go trials');
            ax_rtIncorr = subplot(5,1,3); hold on; title(ax_rtIncorr,'RT on incorrect go trials');
            ax_perf = subplot(5,1,4); hold on; title(ax_perf,'Performance');
            ax_png = subplot(5,1,5); hold on; title(ax_png,'proportion of NoGo choices');
            xlabel(ax_png,xlab)
            linkaxes([ax_delays ax_rtCorr ax_rtIncorr ax_perf ax_png],'x');
            xlim([-0.4 0.4]);
            
            %First: discretise laser onset delays. Given that the total
            %dataset is not uniformly sampled, use clustering.
            numBins = 30;
            [bins,edges] = discretize(d.plotDelay,numBins);
            centroids = edges(1:end-1) + 0.5*(edges(2)-edges(1));
            for clu = 1:numBins
                hx = histogram(ax_delays,d.plotDelay(bins==clu),linspace(-0.5,0.5,1000),'Normalization','count');
                hx.EdgeAlpha=0;
            end
            d.plotDelay = centroids(bins)';
            plot(ax_delays,centroids,0,'r.','Markersize',20)
            
            %
            %             if length(delays) > 20
            % %                 [bins,edges] = discretize(obj.data{n}.laserOnset,12);
            % %                 centroids = edges(1:end-1) + 0.5*(edges(2)-edges(1));
            % %
            % %                 histogram(obj.data{n}.laserOnset,100);
            % %
            % %
            %                 numClusters = 15;
            %                 [IDX, centroids] = kmeans(d.plotDelay, numClusters);
            %                 for clu = 1:numClusters
            %                     hx = histogram(ax_delays,d.plotDelay(IDX==clu),linspace(-0.5,0.5,1000),'Normalization','count');
            %                     hx.EdgeAlpha=0;
            %                 end
            %                 plot(ax_delays,centroids,0,'r.','Markersize',30)
            %                 d.plotDelay = centroids(IDX);
            %             else
            %                 centroids = delays;
            %                 histogram(ax_delays,d.plotDelay,100,'Normalization','count');
            %                 title(ax_delays,'Preset delay (not measured)');
            %             end
            %             xlim([min(centroids) max(centroids)] + [-0.05 +0.05]);
            
            %Second, calculate non-laser measures of
            %choices/performance/reaction time. Only looking at trials
            %where CL>CR or CR>CL
            %             d = obj.data{n};
            %             d = getrow(d, (d.stimulus(:,1)>0 | d.stimulus(:,2)>0) & d.stimulus(:,1)~=d.stimulus(:,2) );
            
            nl = struct;
            nl.idx = d.laserType==0;
            nl.reactiontimeCorr = mean( d.RT(d.feedbackType==1 & nl.idx) ); %avg RT on correct go trials
            nl.reactiontimeIncorr = mean( d.RT(d.feedbackType==0 & d.response<3 & nl.idx) ); %avg RT on incorrect go trials
            nl.performance = mean(d.feedbackType(nl.idx)); %Overall performance on nonlaser trials
            nl.pNoGo = mean(d.response(nl.idx)==3);
            
            
            %Now calculate the value of these at each timepoint of delay
            
            %             laserOnsets = unique(d.laserOnset(~isnan(d.laserOnset)));
            laserLocations = [1 5];
            locationLabels = {'contra V1 (LV1 for CR>CL)+(RV1 for CL>CR)','contra M2 (LM2 for CR>CL)+(RM2 for CL>CR)'};
            locationCols = [      0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.9290    0.6940    0.1250;
                0.4940    0.1840    0.5560;
                0.4660    0.6740    0.1880];
            jitter = linspace(-0.01,0.01,length(laserLocations))/10;
            
            %Also add a field for movement time relative to laser onset
            
            
            for laserPos = 1:length(laserLocations)
                %Select contrasts where contra>ipsi for a given
                %inactivation location
                e = getrow(d, (d.areaIdx == laserLocations(laserPos) & d.stimulus(:,2)>d.stimulus(:,1)) | (d.areaIdx == laserLocations(laserPos)+1 & d.stimulus(:,1)>d.stimulus(:,2)) );
                
                las = struct;
                for onset = 1:length(centroids)
                    las.idx = e.plotDelay==centroids(onset);
                    
                    rt = e.RT(e.feedbackType==1 & las.idx);
                    las.reactiontimeCorr(onset,1) = mean( rt );
                    las.reactiontimeCorr_err(onset,1) = std( rt )/sqrt(length(rt));
                    
                    rt = e.RT(e.feedbackType==0 & e.response<3 & las.idx);
                    las.reactiontimeIncorr(onset,1) = mean( rt );
                    las.reactiontimeIncorr_err(onset,1) = std( rt )/sqrt(length(rt));
                    
                    [p,pci]=binofit(sum(e.feedbackType(las.idx)), length(e.feedbackType(las.idx)));
                    las.performance(onset,1) = p;
                    las.performance_ci(onset,:) = pci;
                    
                    [p,pci]=binofit(sum(e.response(las.idx)==3),length(e.response(las.idx)));
                    las.pNoGo(onset,1) = p;
                    las.pNoGo_ci(onset,:) = pci;
                end
                
                %Plot
                h=[];
                h(1) = errorbar(ax_rtCorr,centroids + jitter(laserPos),las.reactiontimeCorr,las.reactiontimeCorr_err);
                h(2) = errorbar(ax_rtIncorr,centroids + jitter(laserPos),las.reactiontimeIncorr,las.reactiontimeIncorr_err);
                h(3) = errorbar(ax_perf,centroids + jitter(laserPos),las.performance,las.performance-las.performance_ci(:,1),las.performance_ci(:,2)-las.performance);
                h(4) = errorbar(ax_png,centroids + jitter(laserPos),las.pNoGo,las.pNoGo-las.pNoGo_ci(:,1),las.pNoGo_ci(:,2)-las.pNoGo);
                set(h,'Color',locationCols(laserPos,:),'Capsize',0,'Marker','.','MarkerSize',20,'LineStyle','-');
            end
            legend(ax_png,locationLabels,'Location','north');
            legend boxoff;
            
            %Plot nonLaser values as lines
            l=[];
            l(1) = line(ax_rtCorr,get(ax_rtCorr,'XLim'),[1 1]*nl.reactiontimeCorr);
            l(2) = line(ax_rtIncorr,get(ax_rtIncorr,'XLim'),[1 1]*nl.reactiontimeIncorr);
            l(3) = line(ax_perf,get(ax_perf,'XLim'),[1 1]*nl.performance);
            l(4) = line(ax_png,get(ax_png,'XLim'),[1 1]*nl.pNoGo);
            set(l,'Color',[1 1 1]*0.5,'LineStyle',':');
            
            %             ylim(ax_rtCorr,[0 0.5]);
            %             ylim(ax_rtIncorr,[0 0.5]);
            
            
        end
        
        function plotEffectOfTime_scatter(obj,n)
            if length(unique(obj.data{n}.laserPower)) > 2
                powers = unique(obj.data{n}.laserPower);
                powers = powers(2:end);
                
                fprintf('<strong> Laser powers available, please select one </strong>\n');
                for i = 1:length(powers)
                    fprintf('%i) %0.1f mW\n',i,powers(i));
                end
                
                p = input('ID=');
                
                if p <= length(powers)
                    obj.data{n} = getrow(obj.data{n}, obj.data{n}.laserPower==0 | obj.data{n}.laserPower==powers(p));
                else
                    error('not valid power option');
                end
            end
            
            
            d = obj.data{n};
            d.RT(d.response==3) = 0;
            
            cont = d.stimulus;
            c_0 = sum(cont,2)==0;
            c_L = cont(:,1) > 0 & cont(:,2)==0;
            c_R = cont(:,2) > 0 & cont(:,1)==0;
            c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
            d.stim = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
            
            figure('color','w','name',obj.names{n});
            laserLocations = {'L V1',1; 'R V1',2; 'L M2', 5;'R M2', 6};
            
            contrastLabels = {'CL only','CR only','CL=CR','C=[0,0]'};
            responseColours = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
            responseLabels = {'Left','Right','NoGo'};
            for laserLoc = 1:size(laserLocations,1)
                for s = 1:4
                    
                    subplot(size(laserLocations,1),4,s + 4*(laserLoc-1));
                    hold on;
                    
                    for r=1:3
                        e = getrow(d,d.areaIdx==laserLocations{laserLoc,2} & d.stim==s & d.response==r);
                        h=plot(e.laserOnset,e.RT,'.','MarkerSize',5);
                        h.Color = responseColours{r};
                    end
                    
                    ezplot('y=x'); xlim([-1 1]*0.4); ylim([0 1]);
                    title({contrastLabels{s},laserLocations{laserLoc,1}});
                    xlabel('laser time'); ylabel('movement time');
                    
                end
            end
            
        end
        
        function plotChoiceEffects(obj,n)
            %Plot model-free delta pL,pR,pNG, and if available compare with model predicted delta pL,pR,pNG
            if length(unique(obj.data{n}.laserPower)) > 2
                powers = unique(obj.data{n}.laserPower);
                powers = powers(2:end);
                
                fprintf('<strong> Laser powers available, please select one </strong>\n');
                for i = 1:length(powers)
                    fprintf('%i) %0.1f mW\n',i,powers(i));
                end
                
                p = input('ID=');
                
                if p <= length(powers)
                    obj.data{n} = getrow(obj.data{n}, obj.data{n}.laserPower==0 | obj.data{n}.laserPower==powers(p));
                else
                    error('not valid power option');
                end
            end
            
            
            dotShape = 'o';
            
            cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
                ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
            
            img=imread('D:\kirkcaldie_brain_BW_outline.png');
            
            %             numSubjects = length(obj.names);
            
            withShuffle = 0;
            %             if nargin>2
            %                 if contains(varargin{:},'withShuffle')
            %                     withShuffle = 1;
            %
            %                 end
            %             end
            %
            fisherExactP = nan(4,6,2);
            
            if withShuffle == 1
                numIter = 1000;
                
                val_stim_diff_null = nan((size(obj.inactivationCoords,1)), 4, 3, numIter);
                
                for iter = 1:numIter
                    shuf = obj.shuffleData(n);
                    
                    cont = shuf.stimulus(:,1:2);
                    c_0 = sum(cont,2)==0;
                    c_L = cont(:,1) > 0 & cont(:,2)==0;
                    c_R = cont(:,2) > 0 & cont(:,1)==0;
                    c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
                    stim = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
                    
                    resp = shuf.response;
                    laser = shuf.laserIdx;
                    resp(stim==0) = [];
                    laser(stim==0) = [];
                    cont(stim==0,:) = [];
                    stim(stim==0)= [];
                    
                    for r = 1:3
                        val = (resp==r);
                        val_stim = nan(size(obj.inactivationCoords,1)+1,4);
                        
                        for site = 1:(size(obj.inactivationCoords,1)+1)
                            for s = 1:4
                                val_stim(site,s) = nanmean(val(laser==(site-1) & stim==s));
                            end
                        end
                        
                        val_stim_diff_null(:,:, r, iter) = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:));
                        
                    end
                    
                end
                
                
                %statistical test
                areaLabels = {'Left V1','Right V1','Left S1','Right S1','Left M2','Right M2'};
                stim_labels = {'CL only','CR only','CL=CR','C=[0 0]'};
                resp_labels = {'Left','Right','NoGo'};
                sigString = {'ns','*','**','***','****'};
                disp(obj.names{n});
                
                d = obj.data{n};
                cont = d.stimulus(:,1:2);
                c_0 = sum(cont,2)==0;
                c_L = cont(:,1) > 0 & cont(:,2)==0;
                c_R = cont(:,2) > 0 & cont(:,1)==0;
                c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
                stim = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
                
                for a = 1:6
                    disp(['    ' areaLabels{a}]);
                    
                    for s = 1:4
                        disp(['    ' '    ' stim_labels{s}]);
                        idxNL = d.areaIdx==999 & stim==s;
                        idxL = d.areaIdx==a & stim==s;
                        
                        for r = 1:3
                            
                            diff = mean(d.response(idxL)==r) - mean(d.response(idxNL)==r);
                            
                            null = squeeze(val_stim_diff_null(:,s,r,:)); null=null(:);
                            pValue = min([mean(diff<null) mean(diff>null)]);
                            
                            numStars = 1 + (pValue < 0.0001) + (pValue < 0.001) + (pValue < 0.01) + (pValue < 0.05);
                            
                            disp(['    ' '    ' '    ' resp_labels{r} ' ' sigString{numStars}]);
                        end
                    end
                    
                end
                
            end
            
            %             keyboard;
            
            
            
            
            d = obj.data{n};
            
            cont = d.stimulus(:,1:2);
            c_0 = sum(cont,2)==0;
            c_L = cont(:,1) > 0 & cont(:,2)==0;
            c_R = cont(:,2) > 0 & cont(:,1)==0;
            c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
            %                 c_L = cont(:,1) == 0.54 & cont(:,2)==0;
            %                 c_R = cont(:,2) == 0.54 & cont(:,1)==0;
            %                 c_LR = ( cont(:,1) == 0.54 & cont(:,2)== 0.54 ) & ~c_0;
            stim = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
            stim_labels = {'CL only','CR only','CL=CR','C=[0 0]'};
            %                 stim_labels = {'High CL','High CR','High CL&CR','C=[0 0]'};
            col_labels = {'Chose Left','Chose Right','NoGo'};
            
            resp = d.response;
            laser = d.laserIdx;
            sess = d.sessionID;
            
            resp(stim==0) = [];
            laser(stim==0) = [];
            cont(stim==0,:) = [];
            sess(stim==0)=[];
            stim(stim==0)= [];
            
            
            
            
            %                 %Print Chi2 test for laser effects (collapsing over
            %                 %stimuli, area, choices)
            %                 areaLabels = {'Left V1','Right V1','Left S1','Right S1','Left M2','Right M2'};
            %                 sigString = {'ns','*','**','***','****'};
            %                 disp(obj.names{n});
            %
            %                 try
            %                     for a = 1:6
            %                         disp(['    ' areaLabels{a}]);
            %                         idx = (d.areaIdx==a | d.areaIdx==999);
            %
            %                         %If left hemisphere, then test Left+NoGo/Right on
            %                         %trials with CR>CL
            %                         if mod(a,2) == 1
            %                             idx = idx & (cont(:,2)>cont(:,1));
            %                             [c,~,p]=crosstab(d.response(idx)==2,d.areaIdx(idx));
            %                             txt = 'L+NG/R vs Lsr/NoLsr';
            %                         else
            %                             %If right hemisphere, test Right+Nogo/Left on trials
            %                             %with CL>CR
            %                             idx = idx & (cont(:,1)>cont(:,2));
            %                             [c,~,p]=crosstab(d.response(idx)==1,d.areaIdx(idx));
            %                             txt = 'R+NG/L vs Lsr/NoLsr';
            %                         end
            %
            %                         [~,p,~] = fishertest(c);
            %                         numStars = 1 + (p < 0.0001) + (p < 0.001) + (p < 0.01) + (p < 0.05);
            %                         disp(['    ' '    ' txt ' ' sigString{numStars} ' ' num2str(p)]);
            %                     end
            %                     disp(' ');
            %                 catch
            %                 end
            
            
            
            
            
            fx1=figure('name',obj.names{n},'color','w');
            
            AX = [];
            %Plot actual data
            val_stim_diffR=[];
            for r = 1:3
                val = (resp==r);
                val_stim = nan(size(obj.inactivationCoords,1)+1,4);
                %                     for session=1:max(sess)
                for site = 1:(size(obj.inactivationCoords,1)+1)
                    for s = 1:4
                        val_stim(site,s) = nanmean(val(laser==(site-1) & stim==s));
                    end
                end
                %                     end
                %                     val_stim = nanmean(val_stim,3);
                
                %                     %Change pNoGo to pGo
                %                     if r==3
                %                         val_stim = 1 - val_stim;
                %                     end
                
                %compute change in metric from nonLaser condition
                val_stim_diff = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:));
                
                if r==2
                    val_stim_diffR = val_stim_diff;
                end
                
                ap = obj.inactivationCoords(:,2);
                ml = obj.inactivationCoords(:,1);
                if obj.bilateral_flag==1
                    ap = [ap; ap];
                    ml = [ml; -ml];
                    val_stim_diff = [val_stim_diff; val_stim_diff];
                end
                
                q = quantile(val_stim_diff(:),[0.025 0.975]);
                q = max(abs(q));
                
                
                for s = 1:4
                    %                         AX(r,s)=subplot(3,4,s + 4*(r-1));
                    
                    if withShuffle == 1
                        
                        pValue = nan(size(obj.inactivationCoords,1) , 1);
                        null = val_stim_diff_null(:,s,r,:);
                        null = null(:);
                        
                        for site = 1:(size(obj.inactivationCoords,1))
                            point = val_stim_diff(site,s);
                            %                             null = squeeze(val_stim_diff_null(site,s,r,:));
                            pValue(site) = min([mean(point<=null) mean(point>=null)]);
                        end
                        
                        dotSize = (pValue<0.001)*50 + (pValue<0.01)*50 + (pValue<0.05)*25 + 1;
                        if obj.bilateral_flag==1
                            dotSize = [dotSize; dotSize];
                        end
                    else
                        dotSize = ones(size(val_stim_diff(:,s)))*125;
                    end
                    
                    
                    
                    
                    AX(s,r)=subplot(4,3,3*(s-1) + r);
                    hold on;
                    imX=imagesc(linspace(-4.5,4.5,1000),linspace(3.75,-5.2,1000),img,'Parent',AX(s,r));
                    set(gca,'ydir','normal');
                    set(imX,'alphadata',0.7);
                    h=scatter(ml,ap,dotSize,val_stim_diff(:,s),dotShape,'filled'); axis equal;  drawnow;
                    h.MarkerEdgeColor=[1 1 1]*0.75;
                    %                         caxis([-1 1]*q);
                    caxis([-1 1]*0.8);
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
            
            try
                
                nonlaserparams = median(obj.fitData.params{n},1);
                
                if obj.fitData.perSession == 1
                    
                    warning('Model predictions not appropriate when doing per session fits as it requires averaging session parameters');
                end
                
                laserparams = obj.fitData.deltaParams{n};
                
                
                %                     [ZL,ZR] = obj.getModel(obj.fitData.nonLaserModel,obj.fitData.biasMode,n);
                %                     [ZL_las,ZR_las] = obj.getModel(obj.fitData.laserModel,obj.fitData.biasMode,n);
                nonLaserP = [];
                LaserP = [];
                dR = [];
                for s = 1:4
                    c = cont(stim==s,:);
                    nonLaserP(:,:,s) = mean(obj.calculatePhat(obj.model.ZL,obj.model.ZR,[],nonlaserparams,c),1);
                    for site = 1:size(obj.inactivationCoords,1)
                        LaserP(site,:,s) = mean(obj.calculatePhat(obj.model.deltaZL,obj.model.deltaZR,nonlaserparams,laserparams(site,:),c),1);
                    end
                    
                    %                         %flip pNG to pG
                    %                         nonLaserP(:,3,s) = 1-nonLaserP(:,3,s);
                    %                         LaserP(:,3,s) = 1-LaserP(:,3,s);
                    
                    dP = bsxfun(@minus,LaserP(:,:,s),nonLaserP(:,:,s));
                    dR = [dR, dP(:,2)];
                    %                                                 figure(fx2); subplot(1,4,s);
                    %                                                 plot(val_stim_diff(:,s),dP(:,3),'o'); xlabel('Actual dNG'); ylabel('Pred dNG');
                    %                                                 set(gca,'box','off');
                    %                                                 figure(fx1);
                    %
                    ap = obj.inactivationCoords(:,2);
                    ml = obj.inactivationCoords(:,1);
                    if obj.bilateral_flag==1
                        ap = [ap; ap];
                        ml = [ml; -ml];
                        dP = [dP; dP];
                    end
                    
                    for r = 1:3
                        imX=imagesc(linspace(-4.5,4.5,1000)+11,linspace(3.75,-5.2,1000),img,'Parent',AX(s,r)); set(gca,'ydir','normal');
                        set(imX,'alphadata',0.7);
                        h=scatter(AX(s,r),ml+11,ap,125,dP(:,r),dotShape,'filled'); axis equal;
                        h.MarkerEdgeColor=[1 1 1]*0.75;
                    end
                end
                %                     keyboard;
                %                     titleText = 'Actual & Pred';
                
                
            catch
                
                %                     titleText = '';
            end
            set(get(gcf,'children'),'box','off');
            colormap(cmap);
            
            %                 fx1=figure('name',obj.names{n},'color','w');
            %                 for s=1:4
            %                     %compute ZL and ZR
            %                     z=[];
            %                     for site = 1:(size(obj.inactivationCoords,1)+1)
            %                         r = resp(laser==(site-1) & stim==s);
            %                         r=sum([r==1 r==2 r==3],1)/length(r);
            %
            %                         z(site,1) = log(r(1)/r(3));
            %                         z(site,2) = log(r(2)/r(3));
            %                     end
            %
            %                     dz = bsxfun(@minus,z(2:end,:),z(1,:));
            %                     for i = 1:2
            %                         subplot(4,2,(s-1)*2 + i);
            %                         hold on;
            %                         imX=imagesc(linspace(-4.5,4.5,1000),linspace(3.75,-5.2,1000),img);
            %                         set(gca,'ydir','normal');
            %                         set(imX,'alphadata',0.7);
            %                         h=scatter(ml,ap,200,dz(:,i),'s','filled'); axis equal;
            %                         h.MarkerEdgeColor=[1 1 1]*0.95;
            %                         %                         caxis([-1 1]*q);
            %                         caxis([-1 1]*0.4);
            %                         set(gca,'xtick','','ytick','','xcolor','w','ycolor','w');
            %                         if s==1
            %                             labels={'delta ZL','delta ZR'};
            %                             title(labels{i});
            %                         end
            %
            %                         if i==1
            %                             set(gca,'ycolor','k'); ylabel(stim_labels{s});
            %                         end
            %
            %                     end
            %
            %
            %                 end
            %                 colormap(cmap);
            %                 try
            %                     figure('name',obj.names{n},'color','w');
            %                     ezplot('y=x'); hold on;
            %
            %                     plot(val_stim_diffR(:),dR(:),'ko');
            %                     xlabel('Empirical \Delta p(R)');
            %                     ylabel('Predicted \Delta p(R)');
            %                     title('');
            %                     set(gca,'box','off');
            %                     xlim([-1 1]); ylim([-1 1]); axis square;
            %                 catch
            %                 end
            %
            %                Now plot as a scatter plot of pL vs pR
            
            %                 figure('name',obj.names{n},'color','w');
            %                 subplot(4,3,2);
            %                 kimg=imread('\\basket.cortexlab.net\home\stuff\kirkcaldie_brain_BW.PNG');
            %                 imX=image(-5:1:5,4:-1:-6,kimg); axis square; set(gca,'ydir','normal','xtick','','ytick','','box','off','xcolor','w','ycolor','w');
            %                 set(imX,'alphadata',0.7); hold on;
            %                 numSites = size(obj.inactivationCoords,1);
            %
            %                 %                 cols = get(gca,'ColorOrder');
            %                 %                 cols = [1 0 0;
            %                 %                         0.5 0 0;
            %                 %                         1 0 1;
            %                 %                         0.5 0 0.5;
            %                 %                         0 1 0;
            %                 %                         0 0.5 0;
            %                 %                         0.5 0.5 0.5];
            %                 %                 cols = cols(obj.identifyArea(obj.inactivationCoords),:);
            %                 %                 cols = [linspace(0,1,numSites)' linspace(1,0,numSites)' zeros(numSites,1) ];
            %                 %                 cols(obj.inactivationCoords(:,2)>0,3) = 1;
            %                 %
            %
            %                 cols = [69, 198, 234;
            %                     65, 140, 202;
            %                     234, 155, 196;
            %                     176, 115, 175;
            %                     182, 216, 150;
            %                     252, 200, 102;
            %                     201, 54, 0;
            %                     243, 237, 73;
            %                     125 125 125]/255;
            %
            %                 if obj.bilateral_flag==0
            %                     areaID = [1 1 2 3 3 2 1 1,...
            %                         1 1 2 3 3 2 1 1,...
            %                         5 4 4 3 3 4 4 5,...
            %                         5 5 5 3 3 5 5 5,...
            %                         5 5 6 7 7 6 5 5,...
            %                         6 6 7 7 6 6,...
            %                         7 7 7 7,...
            %                         8 8,...
            %                         9];
            %
            %                     sym = {'o','o'};
            %
            %                 else
            %                     areaID = [3 2 1 1,...
            %                         3 2 1 1,...
            %                         3 4 4 5,...
            %                         3 5 5 5,...
            %                         7 6 5 5,...
            %                         7 6 6,...
            %                         7 7,...
            %                         8];
            %                     sym = {'o','o'};
            %                 end
            %                 cols = cols(areaID,:);
            %
            %
            %                 plotCols = cols;
            %                 ap = obj.inactivationCoords(:,2);
            %                 ml = obj.inactivationCoords(:,1);
            %                 Lidx = ml<0;
            %                 Ridx = ~Lidx;
            %                 if obj.bilateral_flag==1
            %                     ap = [ap; ap];
            %                     ml = [ml; -ml];
            %                     plotCols = [plotCols; plotCols];
            %                 end
            %
            %                 %                 scatter(ml,ap,200,plotCols,sym{1},'filled');
            %                 h=scatter(ml(Lidx),ap(Lidx),200,plotCols(Lidx,:),sym{1});
            %                 set(h,'linewidth',2);
            %                 scatter(ml(Ridx),ap(Ridx),200,plotCols(Ridx,:),sym{2},'filled');
            %                 ylim([-6 6]);
            %                 %
            %                 val_stim = nan(size(obj.inactivationCoords,1)+1,4,3);
            %                 NL = cell(3,4);
            %                 NLci = cell(3,4);
            %                 for r = 1:3
            %                     val = (resp==r);
            %                     for site = 1:(size(obj.inactivationCoords,1)+1)
            %                         for s = 1:4
            %                             val_stim(site,s,r) = nanmean(val(laser==(site-1) & stim==s));
            %                         end
            %                     end
            %
            %                     %nonlaser p and pci
            %                     for s = 1:4
            %                         nonLasertrials = val(laser==0 & stim==s);
            %                         [NL{r,s},NLci{r,s}] = binofit(sum(nonLasertrials),length(nonLasertrials));
            %                     end
            %                 end
            %
            %                 for s = 1:4
            %                     subplot(4,4,s+4); hold on;
            %                     %                                             scatter(val_stim(2:end,s,2),val_stim(2:end,s,1),50,cols,'filled');
            %                     h=scatter(val_stim(1+find(Lidx),s,2),val_stim(1+find(Lidx),s,1),50,cols(Lidx,:),sym{1});
            %                     set(h,'linewidth',2);
            %                     scatter(val_stim(1+find(Ridx),s,2),val_stim(1+find(Ridx),s,1),50,cols(Ridx,:),sym{2},'filled');
            %
            %                     scatter(val_stim(1,s,2),val_stim(1,s,1),100,[0 0 0],'filled');
            %                     %                     l=line(ones(1,2)*NL{2,s},NLci{1,s}); set(l,'linewidth',4,'color',[0 0 0]);
            %                     %                     l=line(NLci{2,s},ones(1,2)*NL{1,s}); set(l,'linewidth',4,'color',[0 0 0]);
            %
            %                     pL = val_stim(1,s,1);
            %                     pR = val_stim(1,s,2);
            %                     pNG = val_stim(1,s,3);
            %                     zl = log(pL/pNG);
            %                     zr = log(pR/pNG);
            %
            %                     pLn = @(pRn)( 1 - pRn*(1+ exp(-zr)));
            %                     ezplot(pLn);
            %                     pLn = @(pRn)( (1 - pRn)/(1+ exp(-zl)));
            %                     a=ezplot(pLn); a.LineStyle='--';
            %
            %                     xlabel('pR'); ylabel('pL'); xlim([0 1]); ylim([0 1]); axis square;
            %                     title(stim_labels{s});
            %
            %                     try
            %                         LaserP; %should throw error if model wasn't fitted
            %                         subplot(4,4,s+8); hold on;
            %                         h=scatter(LaserP(Lidx,2,s),LaserP(Lidx,1,s),50,cols(Lidx,:),sym{1});
            %                         set(h,'linewidth',2);
            %                         scatter(LaserP(Ridx,2,s),LaserP(Ridx,1,s),50,cols(Ridx,:),sym{2},'filled');
            %                         scatter(nonLaserP(1,2,s),nonLaserP(1,1,s),100,[0 0 0],'filled');
            %                         pL = nonLaserP(1,1,s);
            %                         pR = nonLaserP(1,2,s);
            %                         pNG = nonLaserP(1,3,s);
            %                         zl = log(pL/pNG);
            %                         zr = log(pR/pNG);
            %
            %                         pLn = @(pRn)( 1 - pRn*(1+ exp(-zr)));
            %                         ezplot(pLn);
            %                         pLn = @(pRn)( (1 - pRn)/(1+ exp(-zl)));
            %                         a=ezplot(pLn); a.LineStyle='--';
            %
            %                         xlabel('pR'); ylabel('pL'); xlim([0 1]); ylim([0 1]); axis square;
            %                         title('Prediction');
            %                     catch
            %                         warning('Need to fit model to plot model predictions');
            %                     end
            %                 end
            
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
            
            
            %             %Plot results of significance testing as maps
            % %             keyboard;
            %             figure('color','w');
            %             titles={'(R+NG)/L vs Laser/noLaser','(L+NG)/R vs Laser/noLaser'};
            %             for n = 1:numSubjects
            %
            %                 for i=1:2
            %                     subplot(numSubjects,2,2*n -2 + i);
            %                     imagesc((fisherExactP(:,:,i,n)<0.05) + (fisherExactP(:,:,i,n)<0.01) + (fisherExactP(:,:,i,n)<0.001)); axis square;
            %
            %                     if n<numSubjects
            %                         set(gca,'xtick','','ytick','');
            %                     else
            %                         set(gca,'xtick',1:6,'XTickLabel',areaLabels,'xticklabelrotation',45,'ytick',1:4,'yticklabel',stim_labels);
            %                     end
            %
            %                     if n ==1
            %                         title(titles{i});
            %                     end
            %
            %                     if i == 1
            %                         ylabel(obj.names{n});
            %                     end
            %                 end
            %             end
            %
            %             colormap('gray');
            %
            %             %false alarm plots
            %             %             figure('color','w');
            %             %             subplot(1,2,1);
            %             %             plot(FA_L(1:end-1),CR_pL(1:end-1),'o'); xlabel('False alarm left choices'); ylabel('max pL during CR with inactivation');
            %             %             xlim([0 1]); ylim([0 1]); axis square;
            %             %             subplot(1,2,2);
            %             %             plot(FA_R(1:end-1),CL_pR(1:end-1),'o'); xlabel('False alarm right choices'); ylabel('max pR during CL with inactivation');
            %             %             xlim([0 1]); ylim([0 1]); axis square;
        end
        
        function choiceStats(obj)
            cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
                ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
            
            img=imread('D:\kirkcaldie_brain_BW_outline.png');
            stim_labels = {'CL only','CR only','CL=CR','C=[0 0]'};
            resp_labels = {'Left','Right','NoGo'};
            areaLabels = {'Left V1','Right V1','Left S1','Right S1','Left M2','Right M2'};
            
            
            numSubjects = length(obj.names)-1;
            
            numIter = 1000;
            
            %             figure('color','w');
            %             i = 1;
            %             for c = 1:4
            %                 for r = 1:3
            %                     ax(c,r) = subplot(4,3,i);
            %                     hold on;
            %
            %                     set(gca,'xtick','','ytick','');
            %
            %                     if c < 4
            %                         set(gca,'xcolor','w');
            %                     else
            %                         xlabel(resp_labels{r});
            %                     end
            %
            %                     if r > 1
            %                         set(gca,'ycolor','w');
            %                     else
            %                         ylabel(stim_labels{c});
            %                     end
            %
            %                     i = i+1;
            %                 end
            %             end
            
            figure('color','w');
            for a = 1:6
                ax(a) = subplot(3,2,a); hold on;
                title(areaLabels{a});
                set(gca,'xtick',1:3,'XTickLabel',resp_labels,'ytick',1:4,'YTickLabel',stim_labels,'ydir','reverse');
                xlim([0.5 3.5]); ylim([0.5 4.5]); axis square;
            end
            
            shift = linspace(-0.2,0.2,ceil(sqrt(numSubjects))); [shiftY,shiftX]=meshgrid(shift);
            
            for n = 1:numSubjects
                val_stim_diff_null = nan((size(obj.inactivationCoords,1)), 4, 3, numIter);
                
                for iter = 1:numIter
                    shuf = obj.shuffleData(n);
                    
                    cont = shuf.stimulus(:,1:2);
                    c_0 = sum(cont,2)==0;
                    c_L = cont(:,1) > 0 & cont(:,2)==0;
                    c_R = cont(:,2) > 0 & cont(:,1)==0;
                    c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
                    stim = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
                    
                    resp = shuf.response;
                    laser = shuf.laserIdx;
                    resp(stim==0) = [];
                    laser(stim==0) = [];
                    cont(stim==0,:) = [];
                    stim(stim==0)= [];
                    
                    for r = 1:3
                        val = (resp==r);
                        val_stim = nan(size(obj.inactivationCoords,1)+1,4);
                        
                        for site = 1:(size(obj.inactivationCoords,1)+1)
                            for s = 1:4
                                val_stim(site,s) = nanmean(val(laser==(site-1) & stim==s));
                            end
                        end
                        
                        val_stim_diff_null(:,:, r, iter) = bsxfun(@minus,val_stim(2:end,:),val_stim(1,:));
                        
                    end
                    
                end
                
                %statistical test
                
                sigString = {'ns','*','**','***','****'};
                disp(obj.names{n});
                
                d = obj.data{n};
                cont = d.stimulus(:,1:2);
                c_0 = sum(cont,2)==0;
                c_L = cont(:,1) > 0 & cont(:,2)==0;
                c_R = cont(:,2) > 0 & cont(:,1)==0;
                c_LR = ( cont(:,1) == cont(:,2) ) & ~c_0;
                stim = 4*c_0 + 1*c_L + 2*c_R + 3*c_LR;
                
                for s = 1:4
                    for r = 1:3
                        null = squeeze(val_stim_diff_null(:,s,r,:)); null=null(:);
                        
                        idxNL = d.areaIdx==999 & stim==s;
                        
                        diff = [];
                        dotSize = [];
                        coords = [];
                        
                        for a = 1:6
                            idxL = d.areaIdx==a & stim==s;
                            diff(a,1) = mean(d.response(idxL)==r) - mean(d.response(idxNL)==r);
                            pValue = min([mean(diff(a,1)<null) mean(diff(a,1)>null)]);
                            dotSize(a,1) = (pValue<0.001)*10 + (pValue<0.01)*10 + (pValue<0.05)*10 + 1;
                            
                            coords(a,:) = median(d.laser(idxL,:));
                            
                            %Plot dot onto chart
                            %                             shift = linspace(-0.2,0.2,numSubjects);
                            if pValue < 0.001
                                gValue = 0;
                            elseif pValue < 0.01
                                gValue = 0.3;
                            elseif pValue < 0.05
                                gValue = 0.6;
                            else
                                gValue = 1;
                            end
                            
                            h=plot(ax(a),r+shiftX(n),s+shiftY(n),'o','markersize',9,'linewidth',1); drawnow;
                            h.MarkerEdgeColor=[1 1 1]*0.5;
                            h.MarkerFaceColor=[1 1 1]*gValue;
                        end
                        
                        %                         imX=imagesc(linspace(-4.5,4.5,1000) + (n-1)*10,linspace(3.75,-5.2,1000),img,'Parent',ax(s,r));
                        %                         set(gca,'ydir','normal');
                        %                         set(imX,'alphadata',0.7);
                        %                         h=scatter(ax(s,r),coords(:,1) + (n-1)*10,coords(:,2),dotSize,diff,'o','filled'); drawnow;
                        %                         caxis([-1 1]*q);
                        %                         caxis([-1 1]*0.8);
                        %                         ylim([-5.2 3.75]);
                        %                         plot
                        
                    end
                end
                
                %                 colormap(cmap);
                
                for a = 1:6
                    disp(['    ' areaLabels{a}]);
                    
                    for s = 1:4
                        disp(['    ' '    ' stim_labels{s}]);
                        idxNL = d.areaIdx==999 & stim==s;
                        idxL = d.areaIdx==a & stim==s;
                        
                        for r = 1:3
                            
                            diff = mean(d.response(idxL)==r) - mean(d.response(idxNL)==r);
                            
                            null = squeeze(val_stim_diff_null(:,s,r,:)); null=null(:);
                            pValue = min([mean(diff<null) mean(diff>null)]);
                            
                            dotSize = (pValue<0.001)*50 + (pValue<0.01)*50 + (pValue<0.05)*25 + 1;
                            
                            numStars = 1 + (pValue < 0.0001) + (pValue < 0.001) + (pValue < 0.01) + (pValue < 0.05);
                            
                            disp(['    ' '    ' '    ' resp_labels{r} ' ' sigString{numStars}]);
                        end
                    end
                    
                end
            end
            
            
        end
    end
    
    methods (Access=private)
        
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
            eRefs = dat.listExps(name);
            concatenateNew = 0;
            
            try
                e=load(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' name '_EXPREFS.mat']);
                
                new = setdiff(eRefs, e.eRefs(:,1)); %See if there are any new expRefs not included in the saved file.
                %If there are, check their block file exists and recalculate over those
                
                if ~isempty(new)
                    bFiles = dat.expFilePath(new,'Block','m');
                    
                    bFilesExist = cellfun(@(s)exist(s),bFiles)>0;
                    new(~bFilesExist)=[];
                end
                
                if ~isempty(new)
                    eRefs = new;
                    concatenateNew = 1;
                    error('there are new eRefs to search through');
                else
                    eRefs = e.eRefs;
                end
            catch
                warning('Re-calculating now...');
                
                for b = 1:length(eRefs)
                    
                    try
                        p = load(dat.expFilePath(eRefs{b},'parameters','m'));
                        l = laserGLM(eRefs{b}); %check if any lasermanip problems
                        
                        good = 1;
                        
                        %If no laser trials in the block at all, exclude
                        if max(l.data.laserType)==0
                            good = 0;
                        end
                        
                        
                        
                        switch(expt)
                            case 'sparse_bilateral_2D'
                                if p.parameters.numRepeats(23)~=1500
                                    good=0;
                                end
                                
                                if p.parameters.rewardOnStimulus(2,end)~=2.5
                                    good=0;
                                end
                                
                                block = dat.loadBlock(eRefs{b});
                                if block.numCompletedTrials<150
                                    good=0;
                                end
                                
                                if ~any( min(l.data.stimulus,[],2) >0)
                                    good=0;
                                end
                                
                                if size(l.inactivationSite,1) ~= 26
                                    good=0;
                                end
                                
                                
                            case 'sparse_unilateral_1D'
                                if any( min(l.data.stimulus,[],2) >0)
                                    good=0;
                                end
                                
                                %                 p.parameters.rewardOnStimulus
                                % max(p.parameters.rewardOnStimulus(2,:))
                                %                 if max(p.parameters.rewardOnStimulus(2,:)) ~= 1
                                %                     good=0;
                                %                 end
                                %
                                if length(l.data.response)<150
                                    good=0;
                                end
                                
                                if size(l.inactivationSite,1) < 5
                                    good=0;
                                end
                                
                                
                                
                            case 'sparse_unilateral_2D'
                                if ~any( min(l.data.stimulus,[],2) >0)
                                    good=0;
                                end
                                
                                %                 p.parameters.rewardOnStimulus
                                % max(p.parameters.rewardOnStimulus(2,:))
                                %                 if max(p.parameters.rewardOnStimulus(2,:)) ~= 1
                                %                     good=0;
                                %                 end
                                %
                                if length(l.data.response)<150
                                    good=0;
                                end
                                
                                if size(l.inactivationSite,1) < 45
                                    good=0;
                                end
                                
                                if ~any(l.data.laser(:,2)<0)
                                    good=0;
                                end
                                
                            case 'sparse_uni&bilateral_2D'
                                if ~any( min(l.data.stimulus,[],2) >0)
                                    good=0;
                                end
                                
                                if length(l.data.response)<150
                                    good=0;
                                end
                                
                                if size(l.inactivationSite,1) < 20
                                    good=0;
                                end
                                
                                
                            case 'sparse_unilateral_1D_highRes'
                                if ~any( min(l.data.stimulus,[],2) >0)
                                    good=0;
                                end
                                
                                if length(l.data.response)<150
                                    good=0;
                                end
                                
                                if size(l.inactivationSite,1) < 5
                                    good=0;
                                end
                                
                                if size(l.inactivationSite,1) > 10
                                    good=0;
                                end
                                
                            case 'sparse_unilateral_2D_highRes'
                                
                                good=0;
                                if strcmp(l.sessionLabel,'LineAlongLeftPPC')
                                    good=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightPPC')
                                    good=1;
                                elseif strcmp(l.sessionLabel,'LineAlongLeftVisual')
                                    good=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightVisual')
                                    good=1;
                                elseif strcmp(l.sessionLabel,'LineAlongLeftM2')
                                    good=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightM2')
                                    good=1;
                                end
                                
                            case 'sparse_unilateral_2D_PPC'
                                
                                good=0;
                                if strcmp(l.sessionLabel,'LineAlongLeftPPC')
                                    good=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightPPC')
                                    good=1;
                                end
                                
                            case 'sparse_unilateral_2D_VIS'
                                good=0;
                                
                                if strcmp(l.sessionLabel,'LineAlongLeftVisual')
                                    good=1;
                                elseif strcmp(l.sessionLabel,'LineAlongRightVisual')
                                    good=1;
                                end
                                
                            case 'sparse_unilateral_2D_fullPeriodInactivation'
                                good=0;
                                
                                if length(l.data.response)<150
                                    good=0;
                                end
                                
                                if size(l.inactivationSite,1) < 45
                                    good=0;
                                end
                                
                                if ~any(l.data.laser(:,2)<0)
                                    good=0;
                                end
                                
                                if p.parameters.LaserStimFixedPeriod == 0
                                    good=0;
                                end
                                
                            case 'galvo_unilateral_2D'
                                block = dat.loadBlock(eRefs{b});
                                block.events;
                                
                                if ~any(min(l.data.stimulus,[],2) > 0)
                                    good=0;
                                end
                                
                                if length(l.data.response)<150
                                    good=0;
                                end
                                
                                %Exclude any blocks before 18th june
                                %because of bad dot size/power calibrations
                                if datenum(l.expRef(1:10),'yyyy-mm-dd') < datenum('2017-06-18','yyyy-mm-dd')
                                    good = 0;
                                end
                                
                                if any(l.data.laserOnset ~= 0)
                                    good = 0;
                                end
                                
                                if max(l.data.laserType) > 1
                                    good = 0;
                                end
                                
                                
                            case 'galvo_unilateral_2D_tdelay'
                                block = dat.loadBlock(eRefs{b});
                                block.events;
                                
                                if ~any(min(l.data.stimulus,[],2) > 0)
                                    good=0;
                                end
                                
                                if length(l.data.response)<150
                                    good=0;
                                end
                                
                                if ~any(l.data.laserOnset~=0)
                                    good = 0;
                                end
                                
                                if datenum(l.expRef(1:10),'yyyy-mm-dd') < datenum('2017-07-11','yyyy-mm-dd')
                                    good = 0;
                                end
                                
                                if max(l.data.laserType) > 1
                                    good = 0;
                                end
                                
                                if max(l.data.laserDuration) ~= 1.5
                                    good = 0;
                                end
                                %
                                %                                 if max(l.data.laserPower) ~= 3.25
                                %                                     good = 0;
                                %                                 end
                                %                                 if exist();
                                
                                %                                 if size(l.inactivationSite,1) > 8
                                %                                     good = 0;
                                %                                 end
                                
                                %                                 keyboard;
                            case {'galvo_unilateral_2D_tdelayPulse'}
                                block = dat.loadBlock(eRefs{b});
                                block.events;
                                
                                if ~any(min(l.data.stimulus,[],2) > 0)
                                    good=0;
                                end
                                
                                if length(l.data.response)<150
                                    good=0;
                                end
                                
                                if ~any(l.data.laserOnset~=0)
                                    good = 0;
                                end
                                
                                if datenum(l.expRef(1:10),'yyyy-mm-dd') < datenum('2017-07-11','yyyy-mm-dd')
                                    good = 0;
                                end
                                
                                if max(l.data.laserType) > 1
                                    good = 0;
                                end
                                
                                if max(l.data.laserDuration) ~= 0.025
                                    good = 0;
                                end
                                
                                %laserOnset times should be fairly randomly
                                %distributed
                                if length(unique(l.data.laserOnset)) < 20
                                    good = 0;
                                end
                                
                                if strcmp(l.data,'2017-08-18_1_Nyx')
                                    keyboard
                                end
                                
                                %                                 if p.parameters.laserPower ~= 1.36
                                %                                     good = 0;
                                %                                 end
                                
                            otherwise
                                error('choose something');
                        end
                        
                        if good==1
                            eRefs{b,1}
                        end
                        
                    catch
                        good=0;
                    end
                    
                    eRefs{b,2}=good;
                end
                
                if concatenateNew == 1
                    eRefs = [e.eRefs; eRefs];
                end
                
                save(['\\basket.cortexlab.net\home\omnibus_files\' expt '_' name '_EXPREFS.mat'],'eRefs');
            end
        end
        
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
            p(p<0.0001)=0.0001; %ensure probabilities are never zero
            
            logLik = nanmean(log2(p));
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
                
                if pos(2) <= -2 && abs(pos(1))>1 %vis
                    areaID(t) = (pos(1)<0)*1 + (pos(1)>0)*2;
                elseif pos(2) == 0 && abs(pos(1)) > 1 %s1
                    areaID(t) = (pos(1)<0)*3 + (pos(1)>0)*4;
                elseif 2 <= pos(2) && pos(2) <= 3 %m2
                    areaID(t) = (pos(1)<0)*5 + (pos(1)>0)*6;
                elseif pos(2) <= -2 && abs(pos(1)) < 1 %retrosplenial
                    areaID(t) = (pos(1)<0)*7 + (pos(1)>0)*8;
                else
                    areaID(t) = 9;
                end
                
                if isnan(pos(1)) && isnan(pos(2))
                    areaID(t) = 0;
                end
                
            end
        end
        
    end
end