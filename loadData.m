function D = loadData(expRef)
if ~dat.expExists(expRef)
    error('expRef does not exist');
end

expPath = dat.expPath(expRef,'expInfo', 'm');
analysesFile = fullfile('\\basket.cortexlab.net\homes\peterzh\analyses',[expRef '_analyses.mat']);
%Try loading behavioural data
try
    D = load(analysesFile);
    D.response;
catch
    fprintf('analyses.mat file not found for %s, reloading and saving\n',expRef);
    %LOAD BLOCK DATA
    block = dat.loadBlock(expRef);
    D = struct('stimulus',[],'response',[],'repeatNum',[],'feedbackType',[],'RT',[]);

    if isfield(block,'events') %SIGNALS
        numT = length(block.events.responseValues);
        D.stimulus = [block.events.contrastLeftValues' block.events.contrastRightValues'];
        D.response = block.events.responseValues';
        D.response(D.response==0) = 3;
        D.response(D.response==1) = 2;
        D.response(D.response==-1) = 1;

        D.repeatNum = block.events.repeatNumValues';
        D.feedbackType = block.events.feedbackValues';

        goCueTime = block.events.stimulusOnTimes';
        goCueTime = goCueTime(1:length(D.response));
        respTime = block.events.responseTimes';
        respTime = respTime(1:length(D.response));
        D.RT = respTime - goCueTime;
        
        %load wheel position
        pos = block.inputs.wheelValues';
        t = block.inputs.wheelTimes';
        
    elseif isfield(block,'trial') %CHOICEWORLD
        trials = block.trial;
        numT = block.numCompletedTrials;
        
        for t=1:numT
            D.stimulus(t,:) = trials(t).condition.visCueContrast';
            D.response(t,1) = trials(t).responseMadeID';
            D.repeatNum(t,1) = trials(t).condition.repeatNum;
            D.feedbackType(t,1) = trials(t).feedbackType;
            D.RT(t,1) = trials(t).responseMadeTime-trials(t).interactiveStartedTime;
        end
        
        %load wheel position
        pos = block.inputSensorPositions;
        t = block.inputSensorPositionTimes;
        
        goCueTime = [trials.interactiveStartedTime]';
        respTime = [trials.responseMadeTime]'; goCueTime = goCueTime(1:length(respTime));
    end

    %MEASURE ACTUAL RESPONSE TIMES using nick's code
    Fs = 1000;
    params.makePlots = true;
    [moveOnsets, moveOffsets] = wheel.findWheelMoves3(pos, t, Fs, params);
    
    intStartTime = goCueTime;
    resp = D.response;
    hasTurn = resp==1|resp==2;
    
    moveType = wheel.classifyWheelMoves(t, pos, moveOnsets, moveOffsets, intStartTime(hasTurn), respTime(hasTurn), resp(hasTurn));
    
    clear dm; dm.moveOnsets = moveOnsets; dm.moveOffsets = moveOffsets; dm.moveType = moveType;
    plotWheel(t, pos, dm); drawnow; hold on;
        
    %Go through each response Time and find the nearest onset time
    D.startMoveTime = nan(size(D.response));
    for i = 1:length(D.response)
        if D.response(i)<3
            idx = find((moveOnsets - respTime(i))<=0,1,'last');
            D.startMoveTime(i) = moveOnsets(idx);
            estRT = moveOnsets(idx) - goCueTime(i);
            if abs(D.RT(i) - estRT)>0.5 || estRT < 0
                warning('something is wrong with the reaction time estimate, reverting back to old RT for this trial');
%                 xlim([-1 +1]+D.startMoveTime(i));
%                 plot(respTime(i),0,'r.','markersize',30);
%                 plot(goCueTime(i),0,'g.','markersize',30);
%                 h = input('ACCEPT? [y/n]: ','s');
                
%                 switch(h)
%                     case {'Y','y'}
%                     case {'N','n'}
                        estRT = D.RT(i);
%                     otherwise
%                         error('Not an option')
%                 end
            end
            D.RT(i) = estRT;
        end
    end
    
    %Trim trials
    D = structfun(@(f)(f(1:length(D.response),:)),D,'uni',0);
    
    %Save to analyses folder
    save(analysesFile,'-struct','D');
end
    


%Try loading laser data
try
    D = load(analysesFile);
    D.laser;
%     error('hi');
catch
    fprintf('laserDat file not found for %s, reloading and saving\n',expRef);
    %LOAD BLOCK DATA
    block = dat.loadBlock(expRef);

    if isfield(block,'events') %SIGNALS
        coords = block.events.galvoCoordsValues(:,1:2);
        pos = block.events.galvoPosValues';
        coords = coords(abs(pos),:);
        coords(:,1) = sign(pos).*coords(:,1);
        
        laserType = block.events.laserTypeValues'; 
%         if min(laserType)>0; warn
        coords(laserType==0,:) = nan;
        D.laser = coords;
        D.laserType = laserType;
        D.laserPower = block.events.laserPowerValues';
        
        laserTrials = ~isnan(D.laser(:,1));
        
        D.laserPower(~laserTrials) = 0;
%         try
            D.laserDuration = block.events.laserDurationValues';
            D.laserOnset = block.events.laserOnsetDelayValues';
            
%         catch
%             D.laserDuration = ones(size(D.response))*1.5;
%             D.laserOnset = ones(size(D.response))*0;
%         end
        D.laserDuration(~laserTrials) = 0;

        %Check with galvo log
        gl = load(fullfile(expPath,[expRef '_galvoLog.mat']));
        missingTrials = setdiff(1:length(D.response),gl.trialNum);
        if ~isempty(missingTrials)
            warning('Mismatch in trial count with galvo log');
            keyboard;
            %             
%             D.laser(missingTrials,:) = NaN;
%             D.laserType(missingTrials) = 0;
% 
%             goodTrials = (1:length(D.response))';
%             goodTrials(missingTrials+1) = [];
% 
%             D = structfun(@(f)f(goodTrials,:),D,'uni',0);
%             numT = length(D.response);
        end


        
        %If timeline file exists, measure actual timings of stimulus onset
        %and laser onset
        tl_file = dat.expFilePath(expRef,'Timeline','m');
        if exist(tl_file,'file')
            t = load( tl_file );
            vars = {t.Timeline.hw.inputs.name};
            
            b = struct('block',block);
            b.waterValveTimes = block.outputs.rewardTimes';
            b.waterValveValues = block.outputs.rewardValues';
            b.trialStartTimes = block.events.trialNumTimes';
            b.VisOnTimes = block.events.stimulusOnTimes';
            b.LaserOnIdx = b.waterValveValues == min(b.waterValveValues);
            
            t.timebase = t.Timeline.rawDAQTimestamps';
            t.photodiode = t.Timeline.rawDAQData(:,strcmp(vars,'photoDiode'));
            t.waterValve = t.Timeline.rawDAQData(:,strcmp(vars,'waterValve'));

            %             %Normalise traces
            t.waterValve = t.waterValve-min(t.waterValve);
            t.waterValve = t.waterValve/max(t.waterValve);
            
            %             %             subplot(1,4,[1:2]);
            %             %             plot(t.timebase,t.photodiode,t.timebase,t.waterValve);
            %             %             hold on;
            %             %             xlim([10 15]);
            %
            %             %Find times when the reward valve switched on
            dy = gradient(t.waterValve);
            [~,ix] = findpeaks(dy, 'MinPeakDistance',50, 'MinPeakHeight',0.25, 'NPeaks',length(b.waterValveTimes));
            t.waterValveTimes = t.timebase(ix);
            %             plot(t.waterValveTimes,1,'ro');
            %
            if length(t.waterValveTimes) ~= length(b.waterValveTimes)
                warning('unequal reward time vectors');
                keyboard;
            end
            %
            %Identify times in timeline when the trial started
            %5ms after TTL pulse, measured with oscilloscope (PZH)
            t.laserOnTimes = t.waterValveTimes(b.LaserOnIdx) + 5/1000;
            %
            fig=figure('name',expRef);
            delay = nan(size(D.response));
            for tr = 1:length(t.laserOnTimes)
                trial_idx = t.laserOnTimes(tr)-1 < t.timebase & t.timebase < t.laserOnTimes(tr)+1;
                
                
                time = t.timebase(trial_idx);
                pd = t.photodiode(trial_idx);
                pd_noise = [mean(pd(1:2000)) std(pd(1:2000))];
                pd = (pd - pd_noise(1))/pd_noise(2);
                pd = abs(pd);
                %                 pdy = gradient(pd);
                wt = t.waterValve(trial_idx);
                %                 plot(time,pd,time,wt,t.laserOnTimes(tr),0,'g.','markersize',30);
                
                %                 [~,ix] = findpeaks(abs(pdy), 'MinPeakDistance',50, 'MinPeakHeight',0.4,'NPeaks',1);
                ix = find(pd>10,1,'first');
                t.screenOnTime(tr,1) = time(ix);
                
                delay(tr) = t.laserOnTimes(tr) - t.screenOnTime(tr,1);
                fprintf('%s %d/%d laserDelay %0.2g\n', expRef, tr, length(t.laserOnTimes), delay(tr));
                
                %                 [~,d]=min(abs(delay(tr) - onset_delays));
                %
                if abs(D.laserOnset(tr) - delay(tr)) > 0.2
                    
                    %                     If error is huge, something is wrong, so plot and pause
                    %                                             fprintf('%s %d/%d laserDelay %0.1g -> %0.1g\n', expRef, tr, length(t.laserOnTimes), D.laserOnset(tr),onset_delays(d));
                    
                    plot(time,pd/max(pd),time,wt,time(ix),0,'r.',t.laserOnTimes(tr),0,'g.','markersize',30);
                    title({sprintf('preset delay %0.1g', D.laserOnset(tr)),...
                        sprintf('measured delay %0.5g', delay(tr))});
                    keyboard
                    
                end
                
                D.laserOnset(tr) = delay(tr); %Overwrite laserOnset with actual delay
            end
            hist(delay);
            keyboard;
            fig.delete;
        end
        
    elseif isfield(block,'trial') %CHOICEWORLD
        laserManipFile = dat.expFilePath(expRef, 'laserManip', 'm');
        if exist(laserManipFile,'file') > 0
            L=load(laserManipFile);
            
            if size(L.coordList,1) > 50
                laserType = 1;
            elseif size(L.coordList,1) == 26
                laserType = 2;
            else
                laserType = NaN;
            end
            
            if isfield(L,'coordList_unadjusted')
                for n=1:size(L.laserCoordByTrial,1)
                    if ~isnan(L.laserCoordByTrial(n,1))
                        matchIdx = sum(abs(bsxfun(@minus,L.laserCoordByTrial(n,:),L.coordList)),2)==0;
                        L.laserCoordByTrial(n,:) = [L.coordList_unadjusted(matchIdx,:)];
                    end
                end
            end
            
            D.laser = [L.laserCoordByTrial(:,2) L.laserCoordByTrial(:,1)];
            D.laserType = ones(size(D.response))*laserType;
            D.laserType(isnan(D.laser(:,1))) = 0;
            D.laserPower = 1.5*ones(size(D.response));
            D.laserPower(isnan(D.laser(:,1))) = 0;
        else
            warning('laser data not found');
            D.laser = nan(size(D.stimulus));
            D.laserType = zeros(size(D.response));
            D.laserPower = zeros(size(D.response));
        end 
    end
    
    %Trim
    D = structfun(@(f)(f(1:length(D.response),:)),D,'uni',0);

    %Save to analyses folder    
    save(analysesFile,'-struct','D');
end


end