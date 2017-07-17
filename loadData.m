function D = loadData(expRef)
    block = dat.loadBlock(expRef);
    D = struct;
    
    if isfield(block,'events') %SIGNALS
        type = 'signals';
        
        numT = length(block.events.responseValues);
        D.stimulus = [block.events.contrastLeftValues' block.events.contrastRightValues'];
        D.response = block.events.responseValues'; 
            D.response(D.response==0) = 3;
            D.response(D.response==1) = 2;
            D.response(D.response==-1) = 1;
            
        D.repeatNum = block.events.repeatNumValues';
        D.feedbackType = block.events.feedbackValues';
        
        stimOnTime = block.events.stimulusOnTimes';
        stimOnTime = stimOnTime(1:length(D.response));
        respTime = block.events.responseTimes';
        respTime = respTime(1:length(D.response));
        D.RT = respTime - stimOnTime;
        
        try %Try finding laser data
            coords = block.events.galvoCoordsValues(:,1:2);
            pos = block.events.galvoPosValues';
            coords = coords(abs(pos),:);
            coords(:,1) = sign(pos).*coords(:,1);
            
            laserType = block.events.laserTypeValues';
            coords(laserType==0,:) = nan;
            D.laser = coords;
            D.laserType = laserType;
            D.laserPower = block.events.laserPowerValues';
            
            try
            D.laserOnset = block.events.laserOnsetDelayValues';
            catch
            end
            
            %Load galvo log file. There may be some trials missing in the
            %log, indicating that the galvo machine lagged on those trials.
            %Therefore those missing trials in the log were actually
            %non-laser trials. The trial after that trial should also be
            %removed.
            expPath = dat.expPath(expRef,'expInfo', 'm');
            gl = load(fullfile(expPath,[expRef '_galvoLog.mat']));
            missingTrials = setdiff(1:length(D.response),gl.trialNum);
            if ~isempty(missingTrials)
                D.laser(missingTrials,:) = NaN;
                D.laserType(missingTrials) = 0;
                
                goodTrials = (1:length(D.response))';
                goodTrials(missingTrials+1) = [];
                
                D = structfun(@(f)f(goodTrials,:),D,'uni',0);
                numT = length(D.response);
            end
            disp(['loaded laser data ' expRef]);

        catch
            warning('no laser data');
        end
        
        
    elseif isfield(block,'trial')
        type = 'choiceworld';
        
        trials = block.trial;
        numT = block.numCompletedTrials;
        
        for t=1:numT
            D.stimulus(t,:) = trials(t).condition.visCueContrast';
            D.response(t,1) = trials(t).responseMadeID';
            D.repeatNum(t,1) = trials(t).condition.repeatNum;
            D.feedbackType(t,1) = trials(t).feedbackType;
            D.RT(t,1) = trials(t).responseMadeTime-trials(t).interactiveStartedTime;
        end
        
        try
            L=load(dat.expFilePath(expRef, 'laserManip', 'm'));
            
            try
                L.coordList_unadjusted(1,:); %in bilat expts this field exists
                for n=1:size(L.laserCoordByTrial,1)
                    if ~isnan(L.laserCoordByTrial(n,1))
                        matchIdx = sum(abs(bsxfun(@minus,L.laserCoordByTrial(n,:),L.coordList)),2)==0;
                        L.laserCoordByTrial(n,:) = [L.coordList_unadjusted(matchIdx,:) 0];
                    end
                end
                
            catch
                %                     disp('No coordList_unadjusted detected');
            end
            
            D.laser = [L.laserCoordByTrial(:,2) L.laserCoordByTrial(:,1)];
            D.laserType = nan(size(D.response));
            D.laserPower = nan(size(D.response));
            
            disp(['loaded laser data ' expRef]);
        catch
            warning('no laser data');
        end
    else
        error('unidentified type');
%         keyboard;
    end
    
    D = structfun(@(f)(f(1:numT,:)),D,'uni',0);
    
end