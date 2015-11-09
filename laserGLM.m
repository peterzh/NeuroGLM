classdef laserGLM < GLM
    % Subclass which handles experimental sessions involving optogenetic
    % inactivation. This is required since those sessions have different
    % sessions depending on the location and intensity of the laser
    % inactivation. Therefore trials must be split up accordingly.
    
    properties
    end
    
    methods
        function obj = laserGLM(inputData)
            obj@GLM(inputData);
            
            if isa(inputData,'char')
                L=load(dat.expFilePath(inputData, 'laserManip', 'm'));
                obj.data.laser = L.laserCoordByTrial;
            end
        end
        
        function obj = fit(obj,laser)
            if laser==0
                obj.data = obj.getrow(obj.data,isnan(obj.data.laser(:,1)));
            elseif laser==1
                obj.data = obj.getrow(obj.data,~isnan(obj.data.laser(:,1)));
            end
            obj=fit@GLM(obj);
        end
        
        function obj = fitCV(obj)
            obj.data = obj.getrow(obj.data,isnan(obj.data.laser(:,1)));
            obj=fitCV@GLM(obj);
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
        
    end
end