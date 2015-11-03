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
        function obj = simulateGLM(expRef,trueModelString,trueModelParams)
            obj@GLM(expRef);
            obj = obj.setModel(trueModelString);
            obj.trueModel = trueModelString;
            obj.trueParameters = trueModelParams;
            
%             obj.data.contrast_cond = repmat(obj.data.contrast_cond,1000,1); %multiply dataset
            
            sim_phat = obj.calculatePhat(obj.trueParameters, obj.data.contrast_cond);
            
            choice = @(phat)(sum(rand>[0,cumsum(phat(1:2))]));
            for i = 1:length(sim_phat)
                obj.data.response(i,1) = choice(sim_phat(i,:));
            end
            
            obj.data.repeatNum = ones(length(sim_phat),1);
            obj.expRef = ['Simulation_' trueModelString];
        end
        
    end
    
end