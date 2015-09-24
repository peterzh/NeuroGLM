classdef NeuroGLM
    %NeuroGLM parent class containing the plotting and analysis methods for analying behavioural
    %and neural data using Generalised Linear Models
    %   Detailed explanation goes here
    
    properties (Access=public)
        Block;
    end
    
    methods
        %constructor, loads the data
        function obj = NeuroGLM(expRef)
            obj.Block = dat.loadBlock(expRef);
        end
        
        %plot the behavioural curves
        
        
    end
    
end

