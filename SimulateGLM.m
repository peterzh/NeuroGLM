classdef simulateGLM < GLM
    % Subclass allows easy simulation of models and subsequent fitting.
    % Useful for validating the fitting procedure. For example one can
    % simulate decisions from a model with given parameters, then use the
    % model fitting to determine whether model comparison selects the right
    % model and parameter estimation recovers the true parameters
    
    properties
        trueModel;
        trueParameters;
    end
    
    methods
        function obj = simulateGLM(expRef,trueModelString,trueModelParams,multiplydata)
            obj@GLM(expRef);
            obj = obj.setModel(trueModelString);
            obj.trueModel = trueModelString;
            obj.trueParameters = trueModelParams;
            
            obj.data.contrast_cond = repmat(obj.data.contrast_cond,multiplydata,1);
            
            sim_phat = obj.calculatePhat(obj.trueParameters, obj.data.contrast_cond);
            
            choice = @(phat)(sum(rand>[0,cumsum(phat(1:2))]));
            for i = 1:length(sim_phat)
                obj.data.response(i,1) = choice(sim_phat(i,:));
            end
            
            obj.data.repeatNum = ones(length(sim_phat),1);
        end
        
        function obj = resample(obj)
            sim_phat = obj.calculatePhat(obj.trueParameters, obj.data.contrast_cond);
            choice = @(phat)(sum(rand>[0,cumsum(phat(1:2))]));
            for i = 1:length(sim_phat)
                obj.data.response(i,1) = choice(sim_phat(i,:));
            end
        end
        
        function obj = fitCV(obj,folds)
            %overload fitCV method to specify which of the simulated trials
            %are tested in the leave1out crossvalidation. This is necessary
            %since simulations may involve millions of trials and it would
            %be infeasible to have millions of cv folds
            
            %Crossvalidated fitting
            
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',10000,'MaxIter',2000);
            
            C = cvpartition(length(obj.data.response),'LeaveOut');
            
            if exist('folds','var')
                NUMFOLDS = folds;
            else
                NUMFOLDS = C.NumTestSets;
            end
            
            obj.parameterFits = nan(NUMFOLDS,length(obj.parameterLabels));
            for f=1:NUMFOLDS
                disp(['Model: ' obj.modelString '. Fold: ' num2str(f) '/' num2str(NUMFOLDS)]);
                trainIdx = find(C.training(f)==1);
                testIdx = find(C.test(f)==1);
                
                trainContrasts = obj.data.contrast_cond(trainIdx,:);
                trainResponses = obj.data.response(trainIdx);
                testContrast = obj.data.contrast_cond(testIdx,:);
                testResponse = obj.data.response(testIdx);
                
                [obj.parameterFits(f,:),~,exitflag] = fmincon(@(b) obj.calculateLogLik(b, trainContrasts, trainResponses), obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
                
                if ~any(exitflag == [1,2])
                    obj.parameterFits(f,:) = nan(1,length(obj.parameterLabels));
                end
                
                phat = obj.calculatePhat(obj.parameterFits(f,:), testContrast);
                
                obj.p_hat(f,1)=phat(testResponse);
            end
        end
        
    end
    
end