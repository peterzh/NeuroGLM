classdef Q
    properties (Access=public)
        data;
%         wt;
%         params=struct;
        N=1;
    end
    
       
    methods
        function obj = Q(expRef)
            obj.data=struct;
            
            %             for b = 1:length(expRef)
            %                 block = dat.loadBlock(expRef{b});
            block = dat.loadBlock(expRef);
            trials = block.trial;
            d = struct;
            
            for t=1:block.numCompletedTrials
                d.stimulus(t,:) = trials(t).condition.visCueContrast';
                d.action(t,1) = trials(t).responseMadeID';
                d.repeatNum(t,1) = trials(t).condition.repeatNum;
                d.reward(t,1) = trials(t).feedbackType==1;
                
                try
                    d.DA(t,1) = trials(t).condition.rewardVolume(2)>0 & trials(t).feedbackType==1;
                catch
                    d.DA(t,1)=nan;
                end
            end
            obj.data = d;
            
            if max(obj.data.action) == 3
                error('Not implemented for 3 choice task');
            end
            
%             badIdx = ~isnan(obj.data.DA);

        end
        
        function obj = fitAlpha(obj)
%             options = optiset('iterfun',@optiplotfval);%'solver','NLOPT',);
%             Opt = opti('fun',@obj.objective,'bounds',[0.0001 0.0001],[1 +inf],'x0',[0.9 0.9])%,'options',options);
%             [p_est,~,exitflag,~] = Opt.solve;
            p_est = fmincon(@obj.objective,[0.1 1],[],[],[],[],[0 0],[1 10]);
            
            alpha = p_est(1);
            beta = p_est(2);
            %             obj.params.beta = p_est(2);
            %             obj.params.beta = 1;
            
            obj.plot(alpha,beta);

        end
        
        function plot(obj,alpha,beta)
            p.alpha = alpha;
            p.beta = beta;
            
            w_init = obj.fitWINIT(p.alpha,p.beta);
            p.sL_init = w_init(1);
            p.bL_init = w_init(2);
            p.sR_init = w_init(3);
            p.bR_init = w_init(4);
            
            w = obj.calculateWT(p);
            ph = obj.calculatePHAT(w);
            bpt = obj.calculateLOGLIK(ph)/length(obj.data.action);
%             
            figure;
            h(1)=subplot(3,1,1);
            plot([w{1};w{2}]');
            legend({'sL','bL','sR','bR'});

            title(['\alpha: ' num2str(p.alpha) '\beta: ' num2str(p.beta)])
%             
            xt = obj.x;
%             y = bsxfun(@times,[1 -1 1 -1]',x);
            numTrials = length(obj.data.action);

            for n = 1:numTrials
                QL(n) = w{1}(:,n)'*xt{1}(:,n);
                QR(n) = w{2}(:,n)'*xt{2}(:,n);
            end

            h(2)=subplot(3,1,2);
            ax=plotyy(1:numTrials,QR-QL,1:numTrials,obj.data.action==2);
            set(get(ax(2),'children'),'linestyle','none','marker','.','markersize',30)
            ylabel('QR - QL'); 
            line([0 numTrials],[ 0 0]);
            ax(2).YTickLabel={'L','R'}; ax(2).YTick=[0 1];
            xlabel('Trial number');
            title(['Bits per trial ' num2str(bpt)]);
%             title(['Bits per trial on training set: ' num2str(nLogLik/length(obj.data.action))])
%             
            linkaxes(h,'x');
            subplot(3,1,3);
            plot(QL(obj.data.action==1),QR(obj.data.action==1),'bo',...
                QL(obj.data.action==2),QR(obj.data.action==2),'ro')
%             scatter(QL,QR,[],obj.data.action,'filled'); colorbar; caxis([1 2]);
            legend({'Chose L','Chose R'});
            axis square;
            hold on;
            ezplot('y=x',[min(QL) max(QL)]);
            xlabel('QL');ylabel('QR');
            hold off;
        end
    end
    
    methods (Access=public)  
        function nloglik = objective(obj,p_vec)
            p.alpha = p_vec(1);
            p.beta = p_vec(2);
            
            w_init = obj.fitWINIT(p.alpha,p.beta);
            p.sL_init = w_init(1);
            p.bL_init = w_init(2);
            p.sR_init = w_init(3);
            p.bR_init = w_init(4);
            
            w = obj.calculateWT(p);
            ph = obj.calculatePHAT(w);
            nloglik = obj.calculateLOGLIK(ph);
        end
%         
        function xt = x(obj)
            numTrials = length(obj.data.action);
            xt = {[obj.data.stimulus(:,1)'.^obj.N; ones(1,numTrials)],
                  [obj.data.stimulus(:,2)'.^obj.N; ones(1,numTrials)]};
        end
        
        function [W_init] = fitWINIT(obj,alpha,beta)
            numTrials = size(obj.data.action,1);
            rt = obj.data.reward*beta;
            at = obj.data.action;
            xt = obj.x;
              
            Vt = {nan(2,numTrials),nan(2,numTrials)};
            Bt = {nan(2,2,numTrials),nan(2,2,numTrials)};
            
            Bt{1}(:,:,1) = eye(2);
            Bt{2}(:,:,1) = eye(2);
            Vt{1}(:,1) = [0;0];
            Vt{2}(:,1) = [0;0];

            offset=nan(numTrials,1);
            regressors=nan(numTrials,4);
            
%             y = [xt{1}(:,1); -xt{2}(:,1)];
%             offset(1) = y'*[Vt{1}(:,1); Vt{2}(:,1)];
            offset(1) = Vt{1}(:,1)'*xt{1}(:,1) - Vt{2}(:,1)'*xt{2}(:,1);
%             regressors(1,:) = y'*[Bt{1}(:,:,1) eye(2); eye(2) Bt{2}(:,:,1) ];
            regressors(1,:) = [xt{1}(:,1)'*Bt{1}(:,:,1) -xt{2}(:,1)'*Bt{2}(:,:,1)];
            
            for t = 2:numTrials
                prev.rt = rt(t-1);
                prev.at = at(t-1);

                if prev.at==1
                    prev.xt = xt{1}(:,t-1);
                    prev.nt = alpha*prev.rt*prev.xt;
                    prev.Mt = eye(2) - alpha*prev.xt*prev.xt';
                    
                    Vt{1}(:,t) = prev.nt + prev.Mt*Vt{1}(:,t-1);
                    Bt{1}(:,:,t) = prev.Mt*Bt{1}(:,:,t-1);
                    
                    Vt{2}(:,t) = Vt{2}(:,t-1);
                    Bt{2}(:,:,t) = Bt{2}(:,:,t-1);
                elseif prev.at==2
                    prev.xt = xt{2}(:,t-1);
                    prev.nt = alpha*prev.rt*prev.xt;
                    prev.Mt = eye(2) - alpha*prev.xt*prev.xt';
                    
                    Vt{2}(:,t) = prev.nt + prev.Mt*Vt{2}(:,t-1);
                    Bt{2}(:,:,t) = prev.Mt*Bt{2}(:,:,t-1);
                    
                    Vt{1}(:,t) = Vt{1}(:,t-1);
                    Bt{1}(:,:,t) = Bt{1}(:,:,t-1);
                end
                
%                 y = [xt{1}(:,t); -xt{2}(:,t)];
%                 offset(t) = y'*[Vt{1}(:,t); Vt{2}(:,t)];
                offset(t) = Vt{1}(:,t)'*xt{1}(:,t) - Vt{2}(:,t)'*xt{2}(:,t);
%                 regressors(t,:) = y'*[Bt{1}(:,:,t) eye(2); eye(2) Bt{2}(:,:,t) ];
                regressors(t,:) = [xt{1}(:,t)'*Bt{1}(:,:,t), -xt{2}(:,t)'*Bt{2}(:,:,t)];


            end
            
            X = regressors;
            Y = obj.data.action;
            Y(Y==2)=0;
            g = fitglm(X,Y,'Distribution','binomial','intercept',false,'offset',offset);
            W_init = table2array(g.Coefficients(:,1));

        end
%         
        function wt = calculateWT(obj,p)
            numTrials = size(obj.data.action,1);
            wt = {nan(2,numTrials),nan(2,numTrials)};
            at = obj.data.action;
            wt{1}(:,1) = [p.sL_init; p.bL_init];
            wt{2}(:,1) = [p.sR_init; p.bR_init];
            
            xt = obj.x;
            
            
            for t = 2:numTrials
                prev.rt = obj.data.reward(t-1)*p.beta;% + obj.data.DA(t-1);
                prev.at = at(t-1);

                if prev.at==1
                    prev.xt = xt{1}(:,t-1);
                    wt{1}(:,t) = wt{1}(:,t-1) + p.alpha*(prev.rt - prev.xt'*wt{1}(:,t-1))*prev.xt;
                    wt{2}(:,t) = wt{2}(:,t-1);
                elseif prev.at==2
                    prev.xt = xt{2}(:,t-1);
                    wt{2}(:,t) = wt{2}(:,t-1) + p.alpha*(prev.rt - prev.xt'*wt{2}(:,t-1))*prev.xt;
                    wt{1}(:,t) = wt{1}(:,t-1);
                end

            end

        end 
%         
        function phat = calculatePHAT(obj,wt)
            xt = obj.x;
            
            numTrials = length(obj.data.action);
            
            for n = 1:numTrials
                QL(n) = wt{1}(:,n)'*xt{1}(:,n);
                QR(n) = wt{2}(:,n)'*xt{2}(:,n);
            end
            
            Z = (QL-QR);

            pL = exp(Z)./(1+exp(Z));
            pR = 1-pL;
            
            phat = [pL; pR];
%             phat(phat==0)=0.0001;
        end
        
        function nloglik = calculateLOGLIK(obj,phat)
            
            if size(phat,1) > size(phat,2)
                phat = phat';
            end
            
            for t = 1:length(obj.data.action)
                r =  obj.data.action(t);
                p(t) = phat(r,t);
            end
            nloglik = -sum(log2(p));
        end
    end
end