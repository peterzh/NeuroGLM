classdef laserGLM < GLM
    properties
    end
    
    methods
        function obj = laserGLM(expRef)
            obj@GLM(expRef);
            L=load(dat.expFilePath(expRef, 'laserManip', 'm'));
            obj.data.laser = L.laserCoordByTrial;
        end
        
        function obj = fit(obj,cv_flag)
            obj.data = obj.getrow(obj.data,isnan(obj.data.laser(:,1)));
            obj=fit@GLM(obj,cv_flag);
        end
        
    end
end