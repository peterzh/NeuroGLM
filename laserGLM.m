classdef laserGLM < GLM
    % Subclass which handles experimental sessions involving optogenetic
    % inactivation. This is required since those sessions have different
    % sessions depending on the location and intensity of the laser
    % inactivation. Therefore trials must be split up accordingly.
    
    properties
        ZINPUT
    end
    
    methods
        function obj = laserGLM(inputData)
            obj@GLM(inputData);
            
            if isa(inputData,'char')
                L=load(dat.expFilePath(inputData, 'laserManip', 'm'));
                obj.data.laser = L.laserCoordByTrial;
            end
        end
        
        function obj = setModel(obj,modelString)
            obj.modelString = modelString;
            obj.parameterFits = [];
            obj.parameterStart = [];
            
            switch(modelString)
                case 'Offset' %Model guesses based on the proportion of responses in the data
                    %used as a baseline to compare other models
                    obj.parameterLabels = {'Offset_L','Offset_R'};
                    obj.parameterBounds = [-inf -inf; +inf +inf];
                    obj.ZL = @(P,CL,CR)(P(1)*ones(length(CL),1));
                    obj.ZR = @(P,CL,CR)(P(2)*ones(length(CR),1));
                case 'fullContrasts'
                    uniqueC = unique(obj.data.contrast_cond,'rows');
                    obj.parameterLabels = models.fullContrast('paramLabels',[],[],[],uniqueC);
                    obj.parameterBounds = models.fullContrast('paramBounds',[],[],[],uniqueC);
                    obj.ZL = @(P,CL,CR)(models.fullContrast('L',P,CL,CR,uniqueC));
                    obj.ZR = @(P,CL,CR)(models.fullContrast('R',P,CL,CR,uniqueC));
                case 'fullContrasts-subset'
                    uniqueC = unique(obj.data.contrast_cond,'rows');
                    obj.parameterLabels = models.fullContrast_subset('paramLabels',[],[],[],uniqueC);
                    obj.parameterBounds = models.fullContrast_subset('paramBounds',[],[],[],uniqueC);
                    obj.ZL = @(P,CL,CR)(models.fullContrast_subset('L',P,CL,CR,uniqueC));
                    obj.ZR = @(P,CL,CR)(models.fullContrast_subset('R',P,CL,CR,uniqueC));

                case 'C^N-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0;
                        +inf +inf +inf +inf +inf];
                    obj.ZINPUT = @(data)([data.contrast_cond(:,1) data.contrast_cond(:,2)]);
                    obj.ZL = @(P,INPUT)(P(1) + P(2).*INPUT(:,1).^P(5));
                    obj.ZR = @(P,INPUT)(P(3) + P(4).*INPUT(:,2).^P(5));
                case 'C^N-subset-laser'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','laser-Offset_L','laser-ScaleL_L','laser-Offset_R','laser-ScaleR_R','N'};
                    obj.parameterBounds = [-inf(1,8) 0;
                        +inf(1,8) inf];
                    %Laser indicator variable included in ZL function
                    obj.ZINPUT = @(data)([data.contrast_cond(:,1) data.contrast_cond(:,2) ~isnan(data.laser(:,1))]);
                    obj.ZL = @(P,INPUT)( (1-INPUT(:,3)).*(P(1) + P(2).*INPUT(:,1).^P(9)) + INPUT(:,3).*(P(5) + P(6).*INPUT(:,1).^P(9)) );
                    obj.ZR = @(P,INPUT)( (1-INPUT(:,3)).*(P(3) + P(4).*INPUT(:,2).^P(9)) + INPUT(:,3).*(P(7) + P(8).*INPUT(:,2).^P(9)) );
                case 'Offset-laser'
                    obj.parameterLabels = {'Offset_L','Offset_R','laser-Offset_L','laser-Offset_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf; +inf +inf +inf +inf];
                    obj.ZINPUT = @(data)(~isnan(data.laser(:,1)));
                    obj.ZL = @(P,INPUT)((1-INPUT).*P(1) + INPUT*P(3));
                    obj.ZR = @(P,INPUT)((1-INPUT).*P(2) + INPUT*P(4));
                otherwise
                    error('Model does not exist');
            end
            
            if isempty(obj.parameterStart)
                obj.parameterStart = zeros(1,length(obj.parameterLabels));
            end
        end
        
        function obj = fit(obj)
            
            %Remove trials with repeats
            obj.data = obj.getrow(obj.data,obj.data.repeatNum==1);
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',10000,'MaxIter',2000,'Display','off');
            
            inputs = obj.ZINPUT(obj.data);
            responses = obj.data.response;
            
            [obj.parameterFits,~,exitflag] = fmincon(@(b) obj.calculateLogLik(b, inputs, responses), obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
            if ~any(exitflag == [1,2])
                obj.parameterFits = nan(1,length(obj.parameterLabels));
            end
            
        end
        
        function obj = fitCV(obj)
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            %Remove trials with repeats
            obj.data = obj.getrow(obj.data,obj.data.repeatNum==1);
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',10000,'MaxIter',2000,'Display','off');
            
            C = cvpartition(length(obj.data.response),'LeaveOut');
            obj.parameterFits = nan(C.NumTestSets,length(obj.parameterLabels));
            for f=1:C.NumTestSets
                disp(['Model: ' obj.modelString '. Fold: ' num2str(f) '/' num2str(C.NumTestSets)]);
                trainIdx = find(C.training(f)==1);
                testIdx = find(C.test(f)==1);
                
                z_inputs = obj.ZINPUT(obj.data);
                trainInputs = z_inputs(trainIdx,:);
                testInputs = z_inputs(testIdx,:);

                trainResponses = obj.data.response(trainIdx);
                testResponse = obj.data.response(testIdx);
                
                [obj.parameterFits(f,:),~,exitflag] = fmincon(@(b) obj.calculateLogLik(b, trainInputs, trainResponses), obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
                
                if ~any(exitflag == [1,2])
                    obj.parameterFits(f,:) = nan(1,length(obj.parameterLabels));
                end
                
                phat = obj.calculatePhat(obj.parameterFits(f,:), testInputs);
                
                obj.p_hat(testIdx,1) = phat(testResponse);
            end
        end
    
        
        function D = calculateLaserEffect(obj)
            query_points = linspace(0,0.02,100);
            
            if isempty(obj.ZL)
                error('Set model with setModel(...) method');
            end
            
            noLaser = obj.fit(0);
            Laser = obj.fit(1);
            
            D = struct;
            D.laser = [0;1];
            D.parameterFits = [noLaser.parameterFits; Laser.parameterFits];
            D.Offsets(1,:) = obj.calculatePhat(D.parameterFits(1,:),[0 0]);
            D.Offsets(2,:) = obj.calculatePhat(D.parameterFits(2,:),[0 0]);
            D.Offsets(:,3) = [];
            
            pc_R = nan(length(query_points),2);
            
            for laser=[1,2]
                for leftright=[1,2]
                    
                    phat = [];
                    for c = 1:length(query_points)
                        contrast = [0 0];
                        contrast(leftright) = query_points(c);
                        phat = [phat; obj.calculatePhat(D.parameterFits(laser,:),contrast)];
                    end
                    
                    phat = phat(:,leftright);
                    
                    if leftright==1
                        phat = -phat;
                    end
                    
                    D.Slopes(laser,leftright) = median(gradient(phat,query_points(2)));
                end
            end
            
            D.laserCoord = Laser.data.laser(1,:);
        end
        
        
        function phat = calculatePhat(obj,testParams,inputs)
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end

            zl = obj.ZL(testParams,inputs);
            zr = obj.ZR(testParams,inputs);
            pL = exp(zl)./(1+exp(zl)+exp(zr));
            pR = exp(zr)./(1+exp(zl)+exp(zr));
            pNG = 1 - pL - pR;
            
            phat = [pL pR pNG];

        end
    end
    
    methods (Access={?GLM})
        function logLik = calculateLogLik(obj,testParams,inputs,responses)
            phat = obj.calculatePhat(testParams, inputs);
            logLik = -sum(log( phat(sub2ind(size(phat), [1:length(responses)]', responses)) ));
        end
        
    end
end