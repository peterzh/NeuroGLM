function [D,meta] = loadData(expRef, force)

disp(['Loading ' expRef])
% dbstop if error

%% Try loading cached file, if not, then continue
sessionCacheFile = fullfile('C:\Users\Peter\Desktop\LocalSessionData',[expRef '_cache.mat']);
try
    load(sessionCacheFile);
    
    D.response;
    if nargin == 1
        return
    else
        warning('Force re-calculation');
    end
catch
    warning('Cached data not found for %s. Calculating now...', expRef);
end

%% Check that servers are accessible
baseDirs = {'\\zserver.cortexlab.net\Data\expInfo';
            '\\zubjects.cortexlab.net\Subjects'};

for i = 1:length(baseDirs)
    if ~exist(baseDirs{i},'dir')
        error('%s not accessible',baseDirs{i});
    end
end


%% Search for the block file
[subj,date,num] = dat.parseExpRef(expRef);
blockFiles = fullfile(baseDirs,subj,datestr(date,'YYYY-mm-DD'),num2str(num),[expRef '_Block.mat']);
blockFilesExist = cellfun(@(f)exist(f,'file'),blockFiles); blockFilesExist = blockFilesExist==2;
idx = find(blockFilesExist);
blockFile = blockFiles{idx(1)};

meta = struct;
meta.blockFile = blockFile;

%If notes file exists, extract and save to metadata
notesFile = fullfile(fileparts(blockFile),'notes.txt');
if exist(notesFile,'file')
    meta.notes = fileread(notesFile);
else
    meta.notes = '';
end

%% Load block file and check data exists
try
    load(blockFile);
    if isempty(block); error('No data within block'); end;
    if isfield(block,'events') && length(block.events.responseValues)<10
        error('No/little data within block');
        %         if isempty(block.events.responseValues) || ~isfield(block.events,'contrastLeftValues') || length(block.events.responseValues)<10
        %             error('Bad block');
        %         end
    elseif isfield(block,'trial') && block.numCompletedTrials< 10;
        error('No/little data within block');
    end
catch
    warning('Block file empty or corrupted: %s', blockFile);
    D = struct('stimulus',[],'response',[],'repeatNum',[],'feedbackType',[]);
    meta = struct('blockFile',blockFile,'notes','empty');
    %     save(sessionCacheFile,'D','meta'); disp('Saved empty!');
    return;
end


%% Extract behavioural data
D = struct;

if isfield(block,'events') %SIGNALS
    
    if contains(block.expDef,{'advancedChoiceWorld.m'}) && ~isfield(block.events,'contrastLeftValues')
        error('Old filetype which did not write separate left and right contrast values');
    elseif contains(block.expDef,{'basicChoiceworld.m'})
        D.stimulus = [block.events.contrastsValues' block.events.contrastsValues'];
        D.stimulus( D.stimulus(:,1)>0 ,1)=0;
        D.stimulus( D.stimulus(:,2)<0 ,2)=0;
        D.stimulus(:,1) = -D.stimulus(:,1);
        D.feedbackType = 2*(block.events.hitValues-0.5);
    else
        D.stimulus = [block.events.contrastLeftValues' block.events.contrastRightValues'];
        D.feedbackType = block.events.feedbackValues';
    end
    
    D.response = block.events.responseValues';
    D.response(D.response==0) = 3;
    D.response(D.response==1) = 2;
    D.response(D.response==-1) = 1;
    
    D.repeatNum = block.events.repeatNumValues';
    
    D.time_stimulusOn = block.events.stimulusOnTimes';
    if isfield(block.events,'goCueTimes')
        D.time_goCue =  block.events.goCueTimes';
    else
        D.time_goCue = D.time_stimulusOn;
    end
    D.time_choiceMade = block.events.responseTimes';
    %
    %             D.RT = respTime - goCueTime;
    %             D.responseTime = respTime;
    %
    %load wheel position
    pos = block.inputs.wheelValues';
    t = block.inputs.wheelTimes';
    wheelGain = block.paramsValues(1).wheelGain;
        
    gl_file = fullfile(fileparts(blockFile),[expRef '_galvoLog.mat']);
    %     tl_file = fullfile(fileparts(blockFile),[expRef '_Timeline.mat']);
    
    if exist(gl_file,'file')
        hasLaserData = 1;
    else
        hasLaserData = 0;
    end
    
elseif isfield(block,'trial') %CHOICEWORLD
    trials = block.trial;
    numT = block.numCompletedTrials;
    
    %             stimulus = arrayfun(@(tr) tr.condition.visCueContrast', trials, 'uni', 0)';
    %             stimulus = cat(1,stimulus{:});
    
    for t=1:numT
        D.stimulus(t,:) = trials(t).condition.visCueContrast';
        D.response(t,1) = trials(t).responseMadeID';
        D.repeatNum(t,1) = trials(t).condition.repeatNum;
        D.feedbackType(t,1) = trials(t).feedbackType;
        D.RT(t,1) = trials(t).responseMadeTime-trials(t).interactiveStartedTime;
        D.responseTime(t,1) = trials(t).responseMadeTime;
    end
    
    %load wheel position
    pos = block.inputSensorPositions;
    t = block.inputSensorPositionTimes;
    
    D.time_stimulusOn = [trials.stimulusCueStartedTime]';
    D.time_goCue = [trials.interactiveStartedTime]';
    D.time_choiceMade = [trials.responseMadeTime]';
    
    laserManipFile = fullfile(fileparts(blockFile),[expRef '_laserManip.mat']);
    if exist(laserManipFile,'file')
        hasLaserData = 1;
    else
        hasLaserData = 0;
    end
    
end

%MEASURE ACTUAL response movement times using nick's code
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
D.time_startMove = nan(size(D.response));
for i = 1:length(D.response)
    if D.response(i) ~= 3
        idx = D.time_goCue(i) < t & t < D.time_goCue(i)+1;
        pos1 = pos(idx);
        
        if ~isempty(pos1)
            t1 = t(idx);
            
            %Transform to common scaling
            pos1 = pos1 - pos1(1);
            pos1 = pos1/max(abs(pos1));
            
            if D.response(i) == 1
                idx = find(pos1 > +0.2,1,'first');
            elseif D.response(i) == 2
                idx = find(pos1 < -0.2,1,'first');
            end
            
            if ~isempty(idx)
                D.time_startMove(i) = t1(idx);
            else
                D.time_startMove(i) = D.time_choiceMade(i);
            end
        else
            D.time_startMove(i) = D.time_choiceMade(i);
        end
        
    end
end

f=figure('color','w', 'name', expRef);

%Wheel position at each stimulus on
wheelBaseline = interp1(t, pos, D.time_stimulusOn);

%Wheel position at relative windows
timeSteps = linspace(-1,1.5,200);
wheelAligned = interp1(t, pos, D.time_stimulusOn + timeSteps);

%subtract wheel position at each event t = 0
wheelAligned = wheelAligned - wheelBaseline;
D.wheel_stimulusOn = wheelAligned;
D.wheel_stimulusOn_timesteps = repmat(timeSteps,size(D.wheel_stimulusOn,1),1);

% %Convert wheel counts to wheel angle
% hd = load(['\\zubjects.cortexlab.net\Code\Rigging\config\' upper(meta.rig) '\hardware.mat']);
% mmFactor = hd.mouseInput.MillimetresFactor;


% keyboard;

figure; axis; hold on;
if any(D.response == 1)
    hx=plot(timeSteps, wheelAligned(D.response == 1,:)); 
    set(hx,'Color',[0 0 1 0.2]);
end
if any(D.response == 2)
    hx=plot(timeSteps, wheelAligned(D.response == 2,:)); 
    set(hx,'Color',[1 0 0 0.2]);
end
if any(D.response == 3)
    hx=plot(timeSteps, wheelAligned(D.response == 3,:)); 
    set(hx,'Color',[0 0 0 0.2]);
end

if any(D.response == 1); 
    plot(timeSteps, nanmean(wheelAligned(D.response == 1,:),1), 'b-', 'linewidth',3); 
end;
if any(D.response == 2); 
    plot(timeSteps, nanmean(wheelAligned(D.response == 2,:),1), 'r-', 'linewidth',3); 
end;

xlabel('Time from stimulus onset'); drawnow;

%Add history terms
for tr = 2:length(D.response)
    D.prev_stimulus(tr,:) = D.stimulus(tr-1,:);
    D.prev_response(tr,1) = D.response(tr-1);
    D.prev_feedback(tr,1) = D.feedbackType(tr-1);
end

%Trim trials
shortest_variable_length = min(structfun(@(f)size(f,1),D));
D = structfun(@(f)(f(1:shortest_variable_length,:)),D,'uni',0);

%Now compute reaction time
D.RT = D.time_startMove - D.time_goCue;

meta.rig = block.rigName;

%If laser data doesnt exist, save and quit
if hasLaserData == 0
    save(sessionCacheFile,'D','meta'); disp('Saved!');
    return;
end

%% Now process laser data
if isfield(block,'events') %SIGNALS
    gl_file = fullfile(fileparts(blockFile),[expRef '_galvoLog.mat']);
    
    coords = block.events.galvoCoordsValues(:,1:2);
    pos = block.events.galvoPosValues';
    coords = coords(abs(pos),:);
    coords(:,1) = sign(pos).*coords(:,1);
    
    laserType = block.events.laserTypeValues';
    %         if min(laserType)>0; warn
    coords(laserType==0,:) = nan;
    D.laserCoord = coords;
    D.laserType = laserType;
    D.laserPower = block.events.laserPowerValues';
    
    %If data was acquired before 16th August 2017, then the laserPower
    %supplied was actually laserVoltage due to lack of calibration.
    %Therefore correct these values manually
    if date < datenum('2017-08-16')
        oldpower = [2 2.1   3   3.25 5];
        newpower = [0.1 0.18  0.5   0.86 1.36];
        D.laserPower = arrayfun(@(e) newpower(e==oldpower), D.laserPower);
    end
    
    laserTrials = ~isnan(D.laserCoord(:,1));
    
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
        %                         keyboard;
        %
        D.laserCoord(missingTrials,:) = NaN;
        D.laserType(missingTrials) = 0;
        
        %         %Remove missing trials, and trial after
        %         goodTrials = (1:length(D.response))';
        %         goodTrials([missingTrials]) = [];
        %
        %         D = structfun(@(f)f(goodTrials,:),D,'uni',0);
        numT = length(D.response);
    end
    %         end
    
    %If timeline file exists, measure actual timings of stimulus onset
    %and laser onset
    tl_file = fullfile(fileparts(blockFile),[expRef '_Timeline.mat']);
    if exist(tl_file,'file')
        t = load( tl_file );
        b = struct('block',block);
        
        vars = {t.Timeline.hw.inputs.name};
        column = [t.Timeline.hw.inputs.arrayColumn];
        
        t.timebase = t.Timeline.rawDAQTimestamps';
        t.photodiode = t.Timeline.rawDAQData(:,column(strcmp(vars,'photoDiode')));
        
        numTrials = length(D.response);
        
        b.trialStartTimes = block.events.trialNumTimes(1:numTrials)';
        b.VisOnTimes = block.events.stimulusOnTimes(1:numTrials)';
        
        if any(contains(vars,'waterValve')) %Old data recorded laser TTL through waterValve channel
            t.TTL = t.Timeline.rawDAQData(:,column(strcmp(vars,'waterValve')));
            ttlvals = block.outputs.rewardValues';
            %                 b.TTLTimes = block.outputs.rewardTimes(ttlvals==min(ttlvals))';
            b.TTL_onsets = block.outputs.rewardTimes(1:numTrials)';
        elseif any(contains(vars,'TTL'))
            t.TTL = t.Timeline.rawDAQData(:,column(strcmp(vars,'TTL')));
            try
                b.TTL_onsets = block.outputs.digitalTTLTimes';
                b.TTL_onsets = b.TTL_onsets(1:2:end);
                b.TTL_onsets = b.TTL_onsets(1:numTrials);
            catch
                b.TTL_onsets = block.outputs.TTLTimes(1:numTrials)';
            end
            b.OnsetDelay = (b.TTL_onsets - b.VisOnTimes);
            b.IntendedOnsetDelay = D.laserOnset;
            b.IntendedOnsetDelay = b.IntendedOnsetDelay(1:numTrials);
            %                                 figure; plot(b.IntendedOnsetDelay,b.OnsetDelay - b.IntendedOnsetDelay,'ro');
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
        
        t.TTL_onsets = t.TTL_onsets(1:numTrials);
        b_t_delay = median(t.TTL_onsets - b.TTL_onsets);
        t.trialStartTimes = b.trialStartTimes + b_t_delay;
        t.expScreenOnTime = b.VisOnTimes(1:length(D.response)) + b_t_delay; %Get expected time of stimulus onset
        t.expTTLTime = b.TTL_onsets(1:length(D.response)) + b_t_delay; %Get expected time of TTL onset
        
        %Now go through each trial, and identify the time of the laser
        %and stimulus onset
        %                 fig=figure('name',expRef);
%         pd_stimTime = nan(length(D.response), 1);
%         ttl_stimTime = nan(length(D.response), 1);
%         
        pd_threshold = 0.15;
        ttl_threshold = 0.1;
        %             f=figure('color','w','units','normalized','outerposition',[1 -0.5 0.5 1.5]);
        %             ax_pdTraces1 = subplot(5,1,1); hold on; xlim([-1 1]*0.2);
        %             ax_pdTraces2 = subplot(5,1,2); hold on; xlim([-1 1]*0.2);
        %             ax_ttlTraces1 = subplot(5,1,3); hold on; xlim([-1 1]*0.2);
        %             ax_ttlTraces2 = subplot(5,1,4); hold on; xlim([-1 1]*0.2);
        %             linkaxes([ax_pdTraces1 ax_pdTraces2 ax_ttlTraces1 ax_ttlTraces2],'x');
        
        %             timeSlices = cell(length(D.response),1);
        %             pdSlices =  cell(length(D.response),1);
        %             ttlSlices =  cell(length(D.response),1);
        
        f=figure('color','w','name',expRef);
        pdTrace=subplot(3,1,1); hold on;
        ttlTrace=subplot(3,1,2); hold on;
        for tr = 1:length(D.response)
            
            %Get photodiode and ttl trace near time when we expect the stimulus to
            %appear
            trial_idx = t.trialStartTimes(tr,1) < t.timebase & t.timebase < t.trialStartTimes(tr,1)+mean(diff(t.trialStartTimes));
            
            %1) Get time around screen update
%             trial_idx = t.expScreenOnTime(tr,1)-0.3 < t.timebase & t.timebase < t.expScreenOnTime(tr,1)+0.3;
            time = t.timebase(trial_idx);
            pd = t.photodiode(trial_idx);
            pd = pd - pd(1);
            pd = abs(pd);
            pd = pd/max(pd);
            ix = find(pd > pd_threshold,1,'first');
            pd_stimTime(tr,1) = time(ix);
            hx=plot(pdTrace, time-time(ix), pd, 'k-');
            hx.Color=[0 0 0 0.2];
%             imagesc(pdTrace,time-time(ix),tr,pd');
            
            
            %2) Get time around TTL
%             trial_idx = t.expTTLTime(tr,1)-0.3 < t.timebase & t.timebase < t.expTTLTime(tr,1)+0.3;
%             time = t.timebase(trial_idx);
            ttl = t.TTL(trial_idx);
            ix = find(ttl > ttl_threshold,1,'first');
            ttl_stimTime(tr,1) = time(ix);
            hx=plot(ttlTrace, time-time(ix), ttl, 'k-');
            hx.Color=[0 0 0 0.2];

%             imagesc(ttlTrace,time-time(ix),tr,ttl');
            
            fprintf('%d/%d\n',tr,length(D.response));
        end
        
        xlim(pdTrace,[-1 1]*0.2);
        xlim(ttlTrace,[-1 1]*0.2);
        
        
        
        %             %Plot discrepency between the measures
        %             subplot(5,3,13);
        
        f1=subplot(3,4,9);
        histogram(1000*(pd_stimTime - t.expScreenOnTime),100); title('Screen delay from expected');
        f2=subplot(3,4,10);
        histogram(1000*(ttl_stimTime - t.expTTLTime),100); title('TTL delay from expected ');
        f3=subplot(3,4,11);
        histogram(1000*(ttl_stimTime - pd_stimTime),100); title('Laser-Stim delay ');
        f4=subplot(3,4,12);
        plot(b.IntendedOnsetDelay, ttl_stimTime - pd_stimTime, 'ko');
        drawnow;
        %Plot some things
        
        %             print(f,fullfile('C:\Users\Peter\Desktop\LocalSessionFigs',expRef),'-dpdf','-bestfit');
        
        
        %Calculate the average discrepency between expected stimulus time,
        %and actual stimulus time. Replace every 'actual' time with
        %whatever is calculated from the median. This is required because
        %the photodiode measure seems unreliable sometimes
        %         t.screenOnTime = t.expScreenOnTime + median( pd_stimTime - t.expScreenOnTime );
        %         t.laserOnTime = t.expTTLTime + median( ttl_stimTime - t.expTTLTime ) + 5/1000;
        t.screenOnTime = pd_stimTime;
        t.laserOnTime = ttl_stimTime + 5/1000;
        D.laserOnset = t.laserOnTime - t.screenOnTime;
        
        %Remove trials where the discrepancy between real and intended onset delay times (on laser trials) is huge
        try
            laserTrials = find(D.laserType>0);
            err = D.laserOnset(laserTrials) - b.IntendedOnsetDelay(laserTrials);
            badTrials = laserTrials(abs(err) > 0.1);
            goodTrials = setdiff(1:length(D.response),badTrials);
            D = structfun(@(f) f(goodTrials,:), D, 'uni', 0);
            fprintf('%d bad trials out of %d laser trials\n',length(badTrials),length(laserTrials));
        catch me
            disp(me);
        end
        
    end
    
elseif isfield(block,'trial') %CHOICEWORLD
    laserManipFile = fullfile(fileparts(blockFile),[expRef '_laserManip.mat']);
    
    if exist(laserManipFile,'file') > 0
        L=load(laserManipFile);
        
        if isfield(L,'coordList')
            if size(L.coordList,1) > 50
                laserType = 1;
            elseif size(L.coordList,1) == 26 || size(L.coordList,1) == 25
                laserType = 2;
            else
                laserType = 1;
            end
        else
            laserType=1;
        end
        
        if ~isfield(L,'laserCoordByTrial')
            warning('Code cant handle these old cases without laserCoordByTrial. Not saving anything...');
            meta.notes = [meta.notes '. Exclude'];
            D.laserCoord = nan(size(D.stimulus));
            D.laserType = zeros(size(D.response));
            D.laserPower = zeros(size(D.response));
            D.laserDuration = zeros(size(D.response));
            D.laserOnset = zeros(size(D.response));
            
            return;
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
        else
            if laserType==2 % old setup used Z position to vary M-L, therefore shift the Z coord to ML coord
                %                 keyboard
                
                %Load standard coord set
                load('\\zubjects.cortexlab.net\Code\Rigging\ExpDefinitions\Peter\26CoordSet.mat');
                coordSet = fliplr(coordSet); %Flip because coordSet is originally ML then AP columns. I want AP then ML;
                
                %Overwrite coordinates
                for tr = 1:length(D.response)
                    thisCoord = L.laserCoordByTrial(tr,:);
                    if ~isnan(thisCoord(1))
                        idx = find(L.coordList(:,1)==thisCoord(1) & L.coordList(:,3)==thisCoord(3));
                        L.laserCoordByTrial(tr,:) = [coordSet(idx,:) 0];
                    end
                end
            end
        end
        
        
        
        D.laserCoord = [L.laserCoordByTrial(:,2) L.laserCoordByTrial(:,1)];
        D.laserType = ones(size(D.response))*laserType;
        D.laserType(isnan(D.laserCoord(:,1))) = 0;
        D.laserPower = 1.5*ones(size(D.response));
        D.laserPower(isnan(D.laserCoord(:,1))) = 0;
        
        D.laserDuration = ones(size(D.response))*1.5;
        D.laserDuration(isnan(D.laserCoord(:,1))) = 0;
        
        D.laserOnset = zeros(size(D.response));
        
    elseif isfield(block.trial(1).condition,'rewardOnStimulus') && length(block.trial(1).condition.rewardOnStimulus) == 2 %Some sessions have no laser manip file but still issue laser through rewardOnStimulus vector
        
        laser = [];
        for tr = 1:length(D.response)
            laser(tr,1) = block.trial(tr).condition.rewardOnStimulus(2);
        end
        laser = laser>0;
        
        D.laserCoord = zeros(length(laser),2);
        D.laserCoord(~laser,:) = nan;
        D.laserType = ones(size(D.response));
        D.laserType(isnan(D.laserCoord(:,1))) = 0;
        D.laserPower = ones(size(D.response))*1.5;
        D.laserPower(isnan(D.laserCoord(:,1))) = 0;
        
        D.laserDuration = ones(size(D.response))*1.5;
        D.laserDuration(isnan(D.laserCoord(:,1))) = 0;
        
        D.laserOnset = zeros(size(D.response));
        
    else
        warning('laser data not found');
        D.laserCoord = nan(size(D.stimulus));
        D.laserType = zeros(size(D.response));
        D.laserPower = zeros(size(D.response));
        D.laserDuration = zeros(size(D.response));
        D.laserOnset = zeros(size(D.response));
    end
end

%Trim
D = structfun(@(f)(f(1:length(D.response),:)),D,'uni',0);

%Add region label
laserRegion = nan(size(D.response));
ML = D.laserCoord(:,1);
AP = D.laserCoord(:,2);
side = D.laserType;

regionLabels = {'LeftVIS','RightVIS',...
    'LeftS1','RightS1',...
    'LeftM2','RightM2',...
    'BilatVIS',...
    'BilatS1',...
    'BilatM2',...
    'LeftFrontOutside','RightFrontOutside',...
    'BilatFrontOutside',...
    'LeftBackOutside','RightBackOutside',...
    'BilatBackOutside',...
    'Other'};
%         laserRegion( isnan(ML) & isnan(AP) ) = 1; %Laser off
laserRegion( side == 1 & ( AP <= -2 ) & ( ML < -1 ) )           = 1; %Left V1
laserRegion( side == 1 & ( AP <= -2 ) & ( ML > 1 ) )            = 2; %Right V1
laserRegion( side == 1 & ( AP == 0 ) & ( ML < 0 ) )             = 3; %Left S1
laserRegion( side == 1 & ( AP == 0 ) & ( ML > 0 ) )             = 4; %Right S1
laserRegion( side == 1 & ( 2 <= AP & AP <= 3 ) & ( ML < 0 ) )   = 5; %Left M2
laserRegion( side == 1 & ( 2 <= AP & AP <= 3 ) & ( ML > 0 ) )   = 6; %Right M2

laserRegion( side == 2 & ( AP <= -2 ) & ( ML > 1 ) )            = 7; %Bilat V1
laserRegion( side == 2 & ( AP == 0 ) & ( ML > 0 ) )             = 8; %Bilat S1
laserRegion( side == 2 & ( 2 <= AP & AP <= 3 ) & ( ML > 0 ) )   = 9; %Bilat M2

%Add control regions
laserRegion( side == 1 & ( AP == 3 ) & ( ML == -3.5 ) )         = 10; %Left Front Outside
laserRegion( side == 1 & ( AP == 3 ) & ( ML == 3.5 ) )          = 11; %Right Front Outside
laserRegion( side == 2 & ( AP == 3 ) & ( ML == 3.5 ) )          = 12; %Bilat Front Outside

laserRegion( side == 1 & ( AP == -7 ) & ( ML == -2.5 ) )        = 13; %Left Back Outside
laserRegion( side == 1 & ( AP == -7 ) & ( ML == 2.5 ) )         = 14; %Right Back Outside
laserRegion( side == 2 & ( AP == -7 ) & ( ML == 2.5 ) )         = 15; %Bilat Back Outside

%Mark all other areas
laserRegion( ~isnan(AP) & isnan(laserRegion) ) = 16; %Other regions

D.laserRegion = categorical(laserRegion,1:length(regionLabels),regionLabels);

%Add history terms
for tr = 2:length(D.response)
    D.prev_laser(tr,:) = D.laserCoord(tr-1,:);
    D.prev_laserType(tr,1) = D.laserType(tr-1);
    D.prev_laserPower(tr,1) = D.laserPower(tr-1);
    D.prev_laserRegion(tr,1) = D.laserRegion(tr-1);
end

%Save to analyses folder
save(sessionCacheFile,'D','meta'); disp('Saved!');

close(f);
return;

end