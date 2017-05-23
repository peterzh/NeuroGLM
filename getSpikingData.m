function varargout = getSpikingData(name,varargin)
switch(name)
    case 'nick'
        dirs={'\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150529\cw_1','Laveran Left PPC 1';
            '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150530\cw_1','Laveran Left Cingulate 1';
            '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150601\cw_1','Laveran Right Cingulate 1';
            '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150602\cw_1','Laveran Left V1 1';
            '\\basket.cortexlab.net\data\nick\M150218_NS1LAV\20150528\cw_1','Laveran Left Cingulate 2';
            '\\basket.cortexlab.net\data\nick\M140528_NS1\20141130\cw_1','M140528_NS1 1';
            '\\basket.cortexlab.net\data\nick\M140528_NS1\20141201\cw_1','M140528_NS1 2';
            '\\basket.cortexlab.net\data\nick\M140528_NS1\20141202',     'M140528_NS1 3'};
        
        id = varargin{1};
        cw_dir=dirs{id,1};
        s=load(fullfile(cw_dir,'spikes.mat'));
        
        spikeTimes = cell(length(s.spikes.shanks.singleUnitInds),1);
        for clu = 1:length(spikeTimes)
            spikeTimes{clu} = s.spikes.shanks.units(s.spikes.shanks.singleUnitInds(clu)).spikeTimes;
        end      
        
        %Load behavioural data
        cwl=load(fullfile(cw_dir,'cwLabels.mat'));
        cwe=load(fullfile(cw_dir,'cwEvents.mat'));
        inclTrials = [cwl.cwLabels.repeatNum==1 & cwl.cwLabels.inclTrials]';
        cwl.cwLabels=rmfield(cwl.cwLabels,'contrastIndVals');
        cwl=structfun(@(a)(a(inclTrials)'),cwl.cwLabels,'uni',0);
        cwe=structfun(@(a)(a(inclTrials)),cwe.cwEvents,'uni',0);
        
        %Resampling NoGo trial timings
        GO_idx = cwl.responseMade<3;
        NG_idx = cwl.responseMade==3;
        numNG = sum(NG_idx);
        
        GO_delays = cwe.responseMade(GO_idx) - cwe.goCue(GO_idx);
        cwe.responseMade(NG_idx) = cwe.goCue(NG_idx) + datasample(GO_delays,numNG);
        
        GO_delays = cwe.responseMoveStartTimes(GO_idx) - cwe.goCue(GO_idx);
        cwe.responseMoveStartTimes(NG_idx) = cwe.goCue(NG_idx) + datasample(GO_delays,numNG);
        
        GO_delays =  cwe.feedbackStarted(GO_idx) - cwe.goCue(GO_idx);
        cwe.feedbackStarted(NG_idx) = cwe.goCue(NG_idx) + datasample(GO_delays,numNG);
        
        behav = struct;
        behav.contrast_cond = [cwl.contrastLeft cwl.contrastRight];
        behav.response = cwl.responseMade;
        behav.repeatNum = cwl.repeatNum;
        behav.feedbackType = cwl.feedbackType;
        
        behavT = struct;
        behavT.timestamps = {'stimOnset','goCue','responseMoveStart'};
        behavT.stimOnset = cwe.stimulusCueStarted;
        behavT.goCue = cwe.goCue;
        behavT.responseMoveStart = cwe.responseMoveStartTimes;
        
        otherInfo = {dirs{id,2}};
        
    case 'nick neuropixels'
        dirs = {'\\zserver.cortexlab.net\Data\Subjects\Cori\2016-12-14','ephys_M2','ephys_V1'; %1 2
                '\\zserver.cortexlab.net\Data\Subjects\Cori\2016-12-15','ephys_M2','ephys_V1'; %3 4 ...
                '\\zserver.cortexlab.net\Data\Subjects\Cori\2016-12-16','ephys_AM','ephys_RL'; % 5 6
                '\\zserver.cortexlab.net\Data\Subjects\Cori\2016-12-17','ephys_LM','ephys_PM';
                '\\zserver.cortexlab.net\Data\Subjects\Cori\2016-12-18','ephys_AM','ephys_V1';
                '\\zserver.cortexlab.net\Data\Subjects\Radnitz\2017-01-08','ephys_M2','ephys_V1';
                '\\zserver.cortexlab.net\Data\Subjects\Radnitz\2017-01-09','ephys_M2','ephys_V1';
                '\\zserver.cortexlab.net\Data\Subjects\Radnitz\2017-01-10','ephys_AM','ephys_M2';
                '\\zserver.cortexlab.net\Data\Subjects\Radnitz\2017-01-11','ephys_PM','ephys_RSP';
                '\\zserver.cortexlab.net\Data\Subjects\Radnitz\2017-01-12','ephys_M1','ephys_S1';
                '\\zserver.cortexlab.net\Data\Subjects\Radnitz\2017-01-13','ephys_M1','ephys_V1';
                '\\zserver.cortexlab.net\Data\Subjects\Muller\2017-01-07','ephys_M2','ephys_V1';
                '\\zserver.cortexlab.net\Data\Subjects\Muller\2017-01-08','ephys_AM','ephys_M2';
                '\\zserver.cortexlab.net\Data\Subjects\Muller\2017-01-09','ephys_LM','ephys_PM';
                '\\zserver.cortexlab.net\Data\Subjects\Muller\2017-01-10','ephys_LM','ephys_PM';
                '\\zserver.cortexlab.net\Data\Subjects\Muller\2017-01-11','ephys_M1','ephys_S1';
                '\\zserver.cortexlab.net\Data\Subjects\Muller\2017-01-12','ephys_M1','ephys_V1';
                };
        
        id = varargin{1};
        session = dirs(round(id/2),:);
        
        dat = session{1};
        ephys_shank = session{2+(mod(id,2)==0)};

        
        act = load([dat '\activeData.mat']);
        behav = struct;
        behav.contrast_cond = [act.cweA.contrastLeft act.cweA.contrastRight];
        behav.response = act.cweA.choice;
        behav.repeatNum = act.cweA.repNum;
        behav.feedbackType = act.cweA.feedback;

        behavT = struct;
        behavT.stimOn = act.cwtA.stimOn;
        behavT.beeps = act.cwtA.beeps;
        behavT.moveStart = act.cwtA.moveStart;
        behav = getrow(behav,act.cweA.inclTrials==1);
        behavT = getrow(behavT,act.cweA.inclTrials==1);
        behavT.timestamps = {'stimOn','beeps','moveStart'};

        
        ksDir = fullfile(dat,ephys_shank,'sorting');

        spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));
        spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed
        
        ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
        st = double(ss)/spikeStruct.sample_rate;
        
        %APPLY TIME CORRECTION
        if mod(id,2)==0 %Needs to be aligned
            ephys_timingbase = session{1+(mod(id,2)==0)};
            correctionFile = fullfile(dat,'alignments',['correct_' ephys_shank '_to_' ephys_timingbase '.npy']);
            correction = readNPY(correctionFile);
            
            st = applyCorrection(st,correction);
            
            warning(['aligning ' ephys_shank '_to_' ephys_timingbase]);
        end
        
        if exist(fullfile(ksDir, 'spike_clusters.npy'))
            clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
        else
            clu = spikeTemplates;
        end
        
        
        clusterIDs = unique(clu);
        spikeTimes = cell(length(clusterIDs),1);
        for c = 1:length(clusterIDs)
            spikeTimes{c} = st(clu==clusterIDs(c));
        end
        
        parts = strsplit(dat, '\');
        name = parts{5};
        otherInfo = {[name ' ' ephys_shank(7:end) ' ' num2str(id)]};
        

    case 'sylvia sc'
        info = {'2015-02-23_1_M141002_SS026',...
                'LeftSC?',...
                '\\basket.cortexlab.net\data\nick\M141002_SS026\20150203\cw_1';
                };
            
        id = varargin{1};
        cw_dir=info{id,3};
        s=load(fullfile(cw_dir,'spikes.mat'));
        
        spikeTimes = cell(length(s.spikes.shanks.singleUnitInds),1);
        for clu = 1:length(spikeTimes)
            spikeTimes{clu} = s.spikes.shanks.units(s.spikes.shanks.singleUnitInds(clu)).spikeTimes;
        end      
        
        %Load behavioural data
        cwl=load(fullfile(cw_dir,'cwLabels.mat'));
        cwe=load(fullfile(cw_dir,'cwEvents.mat'));
        inclTrials = [cwl.cwLabels.repeatNum==1 & cwl.cwLabels.inclTrials]';
        cwl.cwLabels=rmfield(cwl.cwLabels,'contrastIndVals');
        cwl=structfun(@(a)(a(inclTrials)'),cwl.cwLabels,'uni',0);
        cwe=structfun(@(a)(a(inclTrials)),cwe.cwEvents,'uni',0);
        
        behav = struct;
        behav.contrast_cond = [cwl.contrastLeft cwl.contrastRight];
        behav.response = cwl.responseMade;
        behav.repeatNum = cwl.repeatNum;
        behav.feedbackType = cwl.feedbackType;
        
        behavT = struct;
        behavT.timestamps = {'stimOnset','goCue','responseMoveStart','responseMade'};
        behavT.stimOnset = cwe.stimulusCueStarted;
        behavT.goCue = cwe.goCue;
        behavT.responseMoveStart = cwe.responseMoveStartTimes;
        behavT.responseMade = cwe.responseMade;
        
        otherInfo = {info{id,2}};
        
    case 'mush v1'
        info = {'2014-11-11_2_M141010_CBECB',...
                'Deep Left V1 1',...
                '\\zserver\Data\multichanspikes\M141010_CBECB\20141111\',...
                'M141010_CBECB_s20141111_all';
                
                '2014-12-02_3_M141010_CBECB',...
                'Deep Left V1 2',...
                '\\zserver\Data\multichanspikes\M141010_CBECB\20141202\',...
                'M141010_CBECB_s20141202_all';
                
                '2014-12-08_3_M141010_CBECB',...
                'Deep Left V1 3',...
                '\\zserver\Data\multichanspikes\M141010_CBECB\20141208\',...
                'M141010_CBECB_s20141208_all';
                
                '2014-12-12_6_M141010_CBECB',...
                'Deep Left V1 4',...
                '\\zserver\Data\multichanspikes\M141010_CBECB\20141212\',...
                'M141010_CBECB_s20141212_all';
                
                '2014-12-13_4_M141010_CBECB',...
                'Deep Left V1 5',...
                '\\zserver\Data\multichanspikes\M141010_CBECB\20141213\',...
                'M141010_CBECB_s20141213_all';
                
                '2014-12-15_2_M141010_CBECB',...
                'Deep Left V1 6',...
                '\\zserver\Data\multichanspikes\M141010_CBECB\20141215\',...
                'M141010_CBECB_s20141215_all';
                
                '2014-12-17_4_M141010_CBECB',...
                'Deep Left V1 7',...
                '\\zserver\Data\multichanspikes\M141010_CBECB\20141217\',...
                'M141010_CBECB_s20141217_all';
                
                '2014-12-18_4_M141010_CBECB',...
                'Deep Left V1 8',...
                '\\zserver\Data\multichanspikes\M141010_CBECB\20141218\',...
                'M141010_CBECB_s20141218_all'};
        
        id = varargin{1};

        res = load([horzcat(info{id,3:4}) '.res.1']);
        clu = load([horzcat(info{id,3:4}) '.clu.1']);
        
        %convert res (counts) to time (sec)
        res = res/30000;
        
        %remove irrelevant value for cluster
        clu(1) = [];
        
        %Go through clusters and pull out spike times
        clusters = unique(clu);
        
        %Remove noise & MUA cluster IDs
        clusters(1:2) = [];
        
        spikeTimes = cell(length(clusters),1);
        for c = 1:length(clusters)
            spikeTimes{c,1} = res(clu==clusters(c));
        end
            
        behav = struct;
               
        %Add behavioural data
        g = GLM(info{id,1});
        behav.contrast_cond = g.data.contrast_cond;
        behav.response = g.data.response;
        behav.repeatNum = g.data.repeatNum;
        behav.feedbackType = g.data.feedbackType;

        numTrials = length(behav.response);
        %Add timestamp information
        b = dat.loadBlock(info{id,1});
        behavT = struct;
        behavT.timestamps = {'onsetToneSoundPlayedTime','stimOnset','responseMade'};
        
        behavT.onsetToneSoundPlayedTime = [b.trial(1:numTrials).onsetToneSoundPlayedTime]';
        behavT.stimOnset = [b.trial(1:numTrials).stimulusCueStartedTime]';
        behavT.responseMade = [b.trial(1:numTrials).responseMadeTime]';

        otherInfo = {info{id,2}};
        
%     case 'all'
%         all = cell(1,9);
%         
%         for n = 1:8
%             [sp,behav,o]=getSpikingData('nick laveran',n);
%             all{n} = {sp,behav,o};
%         end

        
end

varargout = {spikeTimes,behav,behavT,otherInfo};
end