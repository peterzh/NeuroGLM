classdef expressGLM < GLM
    properties
    end
    
    methods
        function obj = expressGLM(mouseName)
            expList = dat.listExps(mouseName);
            
            found=0; offset=0;
            while found==0
                block = dat.loadBlock(expList{end-offset});
                if isempty(block) || ~isstruct(block.trial)
                    disp(['expRef: ' expList{end-offset} 'is empty. Trying earlier...']);
                    offset = offset + 1;
                else
                    found=1;
                end
            end
            
            obj@GLM(expList{end-offset});
            obj = obj.setModel('C^N').fit;
            obj.plotFit;
            set(gcf,'units','normalized','outerposition',[0.2 0.2 0.8 0.8])
            
%             figDir = '//zserver/data/behavior/ChoiceWorld';
        end
    end
    
end