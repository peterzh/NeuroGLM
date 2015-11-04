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
        
        function obj = fit(obj)
            obj.data = obj.getrow(obj.data,isnan(obj.data.laser(:,1)));
            obj=fit@GLM;
        end
        
        function obj = fitCV(obj)
            obj.data = obj.getrow(obj.data,isnan(obj.data.laser(:,1)));
            obj=fitCV@GLM;
        end
        
    end
end