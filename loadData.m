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
%     error('hi');
catch
    fprintf('behavioural postproc file not found for %s, reloading and saving\n',expRef);
    %LOAD BLOCK DATA
    
    try
        block = dat.loadBlock(expRef);
        if isempty(block); dataExists=0; else; dataExists=1; end
    catch
        dataExists=0;
        warning(['corrupted block file ' expRef]);
    end
    
%     D = struct('stimulus',[],'response',[],'repeatNum',[],'feedbackType',[],'RT',[]);
    
    if dataExists
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
        
        %Add history terms
        for tr = 2:length(D.response)
            D.prev_stimulus(tr,:) = D.stimulus(tr-1,:);
            D.prev_response(tr,1) = D.response(tr-1);
            D.prev_feedback(tr,1) = D.feedbackType(tr-1);
        end
        
        %MEASURE ACTUAL RESPONSE TIMES using nick's code
        Fs = 1000;
        rawT = t; rawPos = pos; 
        t = rawT(1):1/Fs:rawT(end);
        pos = interp1(rawT, rawPos, t);
%         params.makePlots = false;
%         [moveOnsets, moveOffsets] = wheel.findWheelMoves3(pos, t, Fs, params);
%         
%         intStartTime = goCueTime;
%         resp = D.response;
%         hasTurn = resp==1|resp==2;
%         
%         moveType = wheel.classifyWheelMoves(t, pos, moveOnsets, moveOffsets, intStartTime(hasTurn), respTime(hasTurn), resp(hasTurn));
%         
%         clear dm; dm.moveOnsets = moveOnsets; dm.moveOffsets = moveOffsets; dm.moveType = moveType;
%             plotWheel(t, pos, dm); drawnow; hold on;
%         plot(t,pos);
        
        %Go through each response Time and find the nearest onset time
        D.startMoveTime = nan(size(D.response));
        for i = 1:length(D.response)
            if D.response(i)<3
                idx = goCueTime(i) < t & t < goCueTime(i)+1;
                pos1 = pos(idx);
                t1 = t(idx);
                
                %Transform to common scaling
                pos1 = pos1 - pos1(1);
                pos1 = pos1/max(abs(pos1));
                
                if D.response(i) == 1
                    idx = find(pos1 > +0.1,1,'first');
                elseif D.response(i) == 2
                    idx = find(pos1 < -0.1,1,'first');
                end
                
                if ~isempty(idx)
                    
                    D.startMoveTime(i) = t1(idx);
                    
                    %                 pos1 = abs(pos1);
                    
%                     plot(t1,pos1,'k:.'); title(D.response(i)); hold on
                    
                    %                 xlim([0 +1]+goCueTime(i));
                    %
                    %                 idx = find(moveOnsets < respTime(i),1,'last');
                    % %                 idx = find((moveOnsets - respTime(i))<=0,1,'last');
                    %                 D.startMoveTime(i) = moveOnsets(idx);
                    %                 estRT = moveOnsets(idx) - goCueTime(i);
                    %                 if abs(D.RT(i) - estRT)>0.5 || estRT < 0
                    %                     warning('something is wrong with the reaction time estimate, reverting back to old RT for this trial');
                    %                                     xlim([-1 +1]+D.startMoveTime(i));
%                     plot(respTime(i),0,'r.','markersize',30);
%                     plot(goCueTime(i),0,'g.','markersize',30);
%                     plot(D.startMoveTime(i),0,'ko','markersize',15);
%                     hold off;
                    %                     %                 h = input('ACCEPT? [y/n]: ','s');
                    %
                    %                     %                 switch(h)
                    %                     case {'Y','y'}
                    %                     case {'N','n'}
                    estRT = D.startMoveTime(i) - goCueTime(i);
                    %                     otherwise
                    %                         error('Not an option')
                    %                     %                 end
                    %                 end
                    D.RT(i) = estRT;
                end
%                 D.RT(i) = estRT;
            end
        end
        
        %Trim trials
        D = structfun(@(f)(f(1:length(D.response),:)),D,'uni',0);
    end
    %Save to analyses folder
    save(analysesFile,'-struct','D');
end



%Try loading laser data
try
    D = load(analysesFile);
    D.laser;
%         error('hi');
catch
    fprintf('laser data not found for %s, reloading and saving\n',expRef);
    %LOAD BLOCK DATA
    if isempty(D.response)
        D.laser = nan(size(D.stimulus));
        D.laserType = zeros(size(D.response));
        D.laserDuration = zeros(size(D.response));
        D.laserPower = zeros(size(D.response));
        D.laserOnset = nan(size(D.response));
    else
        block = dat.loadBlock(expRef);
        
        
        if isfield(block,'events') %SIGNALS
            gl_file = fullfile(expPath,[expRef '_galvoLog.mat']);
            
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
            
            %If data was acquired before 16th August 2017, then the laserPower
            %supplied was actually laserVoltage due to lack of calibration.
            %Therefore correct these values manually
            [~,date,~] = dat.parseExpRef(expRef);
            if date < datenum('2017-08-16')
                oldpower = [2 2.1 3.25 5];
                newpower = [0.1 0.18 0.86 1.36];
                D.laserPower = arrayfun(@(e) newpower(e==oldpower), D.laserPower);
            end
            
            laserTrials = ~isnan(D.laser(:,1));
            
            D.laserPower(~laserTrials) = 0;
            
            %Old data didnt have the laserDuration or laserOnset time fields
            try
                D.laserDuration = block.events.laserDurationValues';
                D.laserOnset = block.events.laserOnsetDelayValues';
            catch
                D.laserDuration = ones(size(D.response))*1.5;
                D.laserOnset = ones(size(D.response))*0;
            end
            D.laserDuration(~laserTrials) = 0;
            
            %Check with galvo log
            %         if exist(gl_file,'file')
            gl = load(gl_file);
            missingTrials = setdiff(1:length(D.response),gl.trialNum);
            if ~isempty(missingTrials)
                warning('Mismatch in trial count with galvo log');
                %                 keyboard;
                %
                D.laser(missingTrials,:) = NaN;
                D.laserType(missingTrials) = 0;
                
                goodTrials = (1:length(D.response))';
                goodTrials(missingTrials+1) = [];
                
                D = structfun(@(f)f(goodTrials,:),D,'uni',0);
                numT = length(D.response);
            end
            %         end
            
            
            
            
            %If timeline file exists, measure actual timings of stimulus onset
            %and laser onset
            tl_file = dat.expFilePath(expRef,'Timeline','m');
            if exist(tl_file,'file')
                t = load( tl_file );
                b = struct('block',block);
                
                vars = {t.Timeline.hw.inputs.name};
                
                t.timebase = t.Timeline.rawDAQTimestamps';
                t.photodiode = t.Timeline.rawDAQData(:,strcmp(vars,'photoDiode'));
                b.trialStartTimes = block.events.trialNumTimes';
                b.VisOnTimes = block.events.stimulusOnTimes';
                
                if any(contains(vars,'waterValve')) %Old data recorded laser TTL through waterValve channel
                    t.TTL = t.Timeline.rawDAQData(:,strcmp(vars,'waterValve'));
                    ttlvals = block.outputs.rewardValues';
                    %                 b.TTLTimes = block.outputs.rewardTimes(ttlvals==min(ttlvals))';
                    b.TTL_onsets = block.outputs.rewardTimes';
                elseif any(contains(vars,'TTL'))
                    t.TTL = t.Timeline.rawDAQData(:,strcmp(vars,'TTL'));
                    try
                        b.TTL_onsets = block.outputs.digitalTTLTimes';
                        b.TTL_onsets = b.TTL_onsets(1:2:end);
                    catch
                        b.TTL_onsets = block.outputs.TTLTimes';
                    end
                    b.OnsetDelay = (b.TTL_onsets - b.VisOnTimes);
                    b.IntendedOnsetDelay = D.laserOnset;
                    b.IntendedOnsetDelay = b.IntendedOnsetDelay(1:length(b.OnsetDelay));
                    figure; plot(b.IntendedOnsetDelay,b.OnsetDelay - b.IntendedOnsetDelay,'ro');
                end
                t.TTL = t.TTL-min(t.TTL);
                t.TTL = t.TTL/max(t.TTL);
                t.TTL = round(t.TTL);
                
                t.TTL_onsets = t.timebase(logical([diff(t.TTL)==1; 0]));
                t.TTL_offsets = t.timebase(logical([diff(t.TTL)==-1; 0]));
                %             pulse_duration = t.TTL_offsets - t.TTL_onsets;
                %             pulse_duration = round(pulse_duration,3);
                %
                %             t.laserOnTimes = t.TTL_onsets + 5/1000;
                
                b_t_delay = mean(t.TTL_onsets - b.TTL_onsets);
                t.trialStartTimes = b.trialStartTimes + b_t_delay;
                
                %Now go through each trial, and identify the time of the laser
                %and stimulus onset
                fig=figure('name',expRef);
                delay = nan(size(D.response));
                pd_relStim = nan(length(D.response), 1001);
                dt = t.timebase(2);

                for tr = 1:length(D.response)
                    %Get idx of trial time period
                    trial_idx = t.trialStartTimes(tr)+0.5 < t.timebase & t.timebase < t.trialStartTimes(tr)+5;
                    time = t.timebase(trial_idx);
                    
                    %Get photodiode measure for that time, and normalise it
                    pd = t.photodiode(trial_idx);
                    preStimWindowIdx = 1:(1/dt);
                    
                    pd = ( pd - mean(pd(preStimWindowIdx)) ) /std(pd(preStimWindowIdx));
                    pd = pd/std(pd(preStimWindowIdx));
                    pd = abs(pd);
                    
                    %                 plot(time,pd)
                    
                    %Find time when photodiode changes significantly. That is
                    %the stimulus onset time
                    ix = find(pd > 15,1,'first');
                    t.screenOnTime(tr,1) = time(ix);
                    
                    %Get the TTL measure
                    ttl = t.TTL(trial_idx);
                    %                 plot(time,ttl);
                    
                    %Get pulse time(s)
                    onsets = time(logical([diff(ttl)==1; 0]));
                    offsets = time(logical([diff(ttl)==-1; 0]));
                    pulse_duration = offsets - onsets;
                    
                    %If multiple pulses (happens when waterValve data was
                    %merged with TTL data in old sessions) then select the
                    %narrowest one
                    t.laserOnTime(tr,1) = onsets(pulse_duration==min(pulse_duration)) + 5/1000;
                    
                    delay(tr) = t.laserOnTime(tr) - t.screenOnTime(tr,1);
                    fprintf('%s %d/%d laserDelay %0.2g\n', expRef, tr, length(D.response), delay(tr));
                    
%                     plot(time,pd/max(pd),time,ttl,t.screenOnTime(tr,1),0,'r.',t.laserOnTime(tr),0,'g.','markersize',30);
                    
                    
                    D.laserOnset(tr) = delay(tr); %Overwrite laserOnset with actual delay
                    
                    %Collate pd measures relative to stimulus onset to verify
                    %at the end that all went well
%                     pd_relStim(tr,:) = pd( (ix-500):(ix+500) )';
                end
%                 h(1)=subplot(3,1,1);
%                 hist(block.events.laserOnsetDelayValues',100); title('Instructed delays');
%                 h(2)=subplot(3,2,3);
%                 hist(b.TTL_onsets - b.VisOnTimes,100); title('Delays in block file');
%                 h(3)=subplot(3,2,4);
%                 plot(D.laserOnset(1:1224),b.OnsetDelay - D.laserOnset(1:1224),'ro')
%                 keyboard;
%                 subplot(3,3,1);
%                 imagesc(pd_relStim);
%                 caxis([0 20])
%                 title('Photodiode relative to detected stimulus onset time'); ylabel('Trial');
%                 set(gca,'xtick',''); xlabel('time');
%                 h(4)=subplot(3,1,3);
%                 hist(delay,100); title('Measured delays from Timeline');
%                 
%                 linkaxes(h,'x');
                
%                 subplot(3,3,4);
%                 hist(t.TTL_onsets - b.TTL_onsets,100); 
%                 title('Timeline TTL onset -  ');
%                 subplot(3,3,5);
%                 hist(t.screenOnTime - b.VisOnTimes,100); 
%                 title('Timeline vis stim onset - Sessions vis stim onset');
%                 drawnow;
                %             keyboard;
                %             fig.delete;
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
                
                if isfield(L,'coordList_unadjusted') && any(L.coordList(:,2)~=L.coordList_unadjusted(:,2))
                    %Correct 3D to 2D position for bilateral sessions
                    
                    for n=1:size(L.laserCoordByTrial,1)
                        if ~isnan(L.laserCoordByTrial(n,1))
                            laserIdx = (L.laserCoordByTrial(n,1) == L.coordList(:,1)) & (L.laserCoordByTrial(n,3) == L.coordList(:,3));
                            L.laserCoordByTrial(n,1:2) = L.coordList_unadjusted(laserIdx,:);
                        end
                    end
                    L.laserCoordByTrial(:,3) = [];
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
    end
    %Save to analyses folder
    save(analysesFile,'-struct','D');
end


end