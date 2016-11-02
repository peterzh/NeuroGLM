classdef bGLM
    properties (Access=public)
        metHastIter = 100000;
        expRef;
        modelString;
        parameterLabels;
        parameterFits;
        parameterStart;
        Zinput;
        ZL;
        ZR;
        data;
        p_hat;
        ContrastDimensions;
    end
    
    properties (Access=public)
        guess_bpt;
        lapseFlag=0;
    end
    
    methods
        function obj = bGLM(inputData)
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
                    end
                catch
                    warning('session data empty');
                end
                
                obj.data = D;
            else
                error('GLM:constructorFail', 'Must pass either an expRef or data struct to GLM constructor');
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
            
            switch(modelString)
                case 'Offset' %Model guesses based on the proportion of responses in the data
                    %used as a baseline to compare other models
                    obj.parameterLabels = {'Offset_L','Offset_R'};
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1)*ones(length(in(:,1)),1));
                    obj.ZR = @(P,in)(P(2)*ones(length(in(:,1)),1));
                case 'C-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R'};
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1)   );
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2)  );
                case 'C'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','ScaleR_L','Offset_R','ScaleL_R','ScaleR_R'};
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1) + P(3).*in(:,2)   );
                    obj.ZR = @(P,in)(P(4) + P(5).*in(:,2) + P(6).*in(:,2)  );
                case 'AFC'
                    obj.parameterLabels = {'Offset','ScaleL','ScaleR'};
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2)*in(:,1) + P(3)*in(:,2));
                    obj.ZR = [];
                case 'AFC_diff'
                    obj.parameterLabels = {'Offset','Scale'};
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2)*(in(:,1)-in(:,2)) );
                    obj.ZR = [];
                case 'AFC_diff_alternate'
                    obj.parameterLabels = {'Offset','Scale'};
                    obj.Zinput = @(D)([D.contrast_cond(:,1) D.contrast_cond(:,2)]);
                    obj.ZL = @(P,in)(P(2)*P(1) + P(2)*(in(:,1)-in(:,2)) );
                    obj.ZR = [];
                otherwise
                    error('Model does not exist');
                    
            end
           
            if isempty(obj.parameterStart)
                obj.parameterStart = zeros(1,length(obj.parameterLabels));
            end

        end
        
        function stan(obj,smObj)
%             switch(modelString)
%                 case '2auc-c-subset'
                    
                    %MULTINOMIAL
                    stanCode = {
                        'data {'
                        '    int K;' %num classes (3)
                        '    int N;' %num trials
                        '    int D;' %num regressors
                        '    int y[N];' %response
                        '    vector[D] x[N];' %input (array of D-vectors)
                        '}'
                        'transformed data {'
                        '    row_vector[D] zeros;'
                        '    zeros = rep_row_vector(0, D);'
                        '}'
                        'parameters {'
                        '    matrix[K - 1,D] beta_raw;'
                        '}'
                        'transformed parameters {'
                        '    matrix[K, D] beta;'
                        '    beta = append_row(beta_raw, zeros);'
                        '}'
                        'model {'
                        '    to_vector(beta_raw) ~ normal(0, 50);'
                        '    for (n in 1:N)'
                        '        y[n] ~ categorical_logit(beta * x[n]);'
                        '}'
                        };
                    
                    dat = struct('K',max(obj.data.response),...
                        'N',length(obj.data.response),...
                        'D',3,...
                        'y',obj.data.response,...
                        'x',[ones(length(obj.data.response),1) obj.data.contrast_cond]);
                    fit=stan('model_code',stanCode,'data',dat,'verbose',true);
                    %takes about 15sec

                    
                    %LOGISTIC 2AFC
                    stanCode = {
                        'data {'
                        '    int numTrials;'
                        '    int y[numTrials];'
                        '    real cl[numTrials];'
                        '    real cr[numTrials];'
                        '}'
                        'parameters {'
                        '    real b;'
                        '    real sl;'
                        '    real sr;'
                        '}'
                        'model {'
                        '    for (n in 1:numTrials)'
                        '        y[n] ~ bernoulli_logit( b + sl*cl[n] + sr*cr[n] );'
                        '}'
                        };
                    
                    D = getrow(obj.data,obj.data.response<3);
                    dat = struct('numTrials',length(D.response),...
                        'y',D.response-1,...
                        'cl',D.contrast_cond(:,1),...
                        'cr',D.contrast_cond(:,2));
                    
                    fit = smObj.sampling('data',dat,'verbose',true);
%                     fit = stan('model_code',stanCode,'data',dat,'verbose',true);
%                     while fit.is_running == 1
%                         pause(0.5);
%                     end
                    
                    keyboard;
                    
                    beta=fit.extract('permute',false).beta;
                    beta(1,:) - beta(3,:)
                    beta(2,:) - beta(3,:)
%                     print(fit);
                    
            
        end

        function p = prior(~,w)
            dim = length(w);
            p = mvnpdf(w,zeros(1,dim),eye(dim)*50);
        end
        
        function l = likelihood(obj,w)
            phat = obj.calculatePhat(w, obj.data.contrast_cond);
            p = (obj.data.response==1).*phat(:,1) + (obj.data.response==2).*phat(:,2) + (obj.data.response==3).*phat(:,3);
            l = prod(p);
        end
        
        function ll = loglikelihood(obj,w)
            phat = obj.calculatePhat(w, obj.data.contrast_cond);
            p = (obj.data.response==1).*phat(:,1) + (obj.data.response==2).*phat(:,2) + (obj.data.response==3).*phat(:,3);
            ll = sum(log(p));
        end
        
        function est_p = posterior(obj)
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            dim = length(obj.parameterLabels);
            %Uses metroplis hastings MCMC            
            w_prev = zeros(1,dim);
            est_p = nan(obj.metHastIter,dim);
            
            for iter = 1:obj.metHastIter
                w_next = mvnrnd(w_prev,eye(dim)*10);
                
                
                logalpha = obj.loglikelihood(w_next) + log(obj.prior(w_next)) - obj.loglikelihood(w_prev) - log(obj.prior(w_prev));
%                 alpha = obj.likelihood(w_next)*obj.prior(w_next) / ...
%                     (obj.likelihood(w_prev)*obj.prior(w_prev));
                alpha = exp(logalpha);
                
                if alpha >= 1 || binornd(1,alpha)==1
                    w_prev = w_next;
                end
                
                est_p(iter,:) = w_prev;
            end
            
            est_p(1:round(obj.metHastIter/10),:) = []; %remove first 10% burn-in
        end
        
        function phat = calculatePhat(obj,testParams,inputs)
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end

            if isempty(obj.ZR) %if a AFC task then no ZR is defined, only pL vs pR
                zl = obj.ZL(testParams,inputs);
                pL = exp(zl)./(1+exp(zl));
                pR = 1 - pL;
                N = length(pL);
                phat = [pL pR zeros(N,1)];
            else 
                zl = obj.ZL(testParams,inputs);
                zr = obj.ZR(testParams,inputs);
                pL = exp(zl)./(1+exp(zl)+exp(zr));
                pR = exp(zr)./(1+exp(zl)+exp(zr));
                pNG = 1 - (pL + pR);
                
                phat = [pL pR pNG];
            end

        end
        
        function post = fit(obj)
            post = obj.posterior;
            figure;
            gplotmatrix(post,[],[],[],[],[],[],[],obj.parameterLabels);
            
            figure;
            imagesc(corrcoef(post)); caxis([-1 1]);
            cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
            colormap(cmap); colorbar;
            dim = size(post,2);
            set(gca,'xtick',1:dim,'ytick',1:dim,'xticklabels',obj.parameterLabels,'yticklabels',obj.parameterLabels);
        end
        
    end
    
end