function varargout = getWidefieldData(name,varargin)
switch(name)
    
    case 'nick widefield'
        %Widefield data is in Timeline time
        %Therefore I should grab the behav data from Timeline directly as
        %well
        
        NUM_COMPONENTS = 10;
        names = {'Moniz','Cori'};
        expRefs = dat.listExps(names);
        id = varargin{1};
        counter = 0;
        
        for n = 1:length(names)
            
            for s = 1:length(expRefs{n})
                
                try
                    g = GLM(expRefs{n}{s});
                    [~,date,session]=dat.parseExpRef(expRefs{n}{s});
                    date = datestr(date,29);
                    
                    %Check if widefield data exists
                    dir = fullfile('\\zserver.cortexlab.net\Data\Subjects\',names{n},date,num2str(session));
                    if ~isempty(g.data.response) && exist(dir) == 7 && mean(g.data.feedbackType==1) > 0.6
                        
                        counter = counter+1;
                        
                        if counter == id
                            
                            disp(expRefs{n}{s});
                            disp(dir);
                            
                            
                            %Load widefield
                            [U,V,t]=quickLoadUVt(dir,NUM_COMPONENTS);
                            W = struct;
                            W.U = U;
                            W.V = V;
                            W.t = t;
                            
                            behav = g.data;
                            
                            %Get timeline
                            tl = load(dat.expFilePath(expRefs{n}{s},'Timeline','master'));
                            b = dat.loadBlock(expRefs{n}{s});
                            
                            t = tl.Timeline.rawDAQTimestamps;
                            t_rewardDelivered = tl.Timeline.rawDAQData(:,find(strcmp({tl.Timeline.hw.inputs.name},'rewardEcho'))); %Reward times in Timeline timebase
                            t_rewardDelivered = t([diff(t_rewardDelivered)>4])';
                            
                            
                            b_stimOnset=nan(b.numCompletedTrials,1);
                            b_goCue=nan(b.numCompletedTrials,1);
                            b_responseMade=nan(b.numCompletedTrials,1);
                            for trial = 1:b.numCompletedTrials
                                b_stimOnset(trial,1)= b.trial(trial).stimulusCueStartedTime;
                                b_goCue(trial,1) = b.trial(trial).onsetToneSoundPlayedTime(1);
                                b_responseMade(trial,1) = b.trial(trial).responseMadeTime;
                            end
                            
                            b_rewardDelivered = b.rewardDeliveryTimes';
                            idx = min([length(b_rewardDelivered) length(t_rewardDelivered)]);
                            figure;
                            while length(b_rewardDelivered) ~= length(t_rewardDelivered)
                                warning('unequal reward timestamps, find missing trial');
                                idx = min([length(b_rewardDelivered) length(t_rewardDelivered)]);
                                badIdx = find(abs(diff(b_rewardDelivered(1:idx)-t_rewardDelivered(1:idx)))>0.1,1,'first')+1;
                                
                                
                                
                                aa=subplot(3,1,1);
                                plot(b_rewardDelivered(1:idx)-t_rewardDelivered(1:idx));
                                
                                
                                if length(b_rewardDelivered) > length(t_rewardDelivered)
                                    b_rewardDelivered(badIdx) = [];
                                else
                                    t_rewardDelivered(badIdx) = [];
                                end
                                
                            end
                            %
                            ab=subplot(3,1,2);
                            plot(b_rewardDelivered-t_rewardDelivered);
                            %                         set(ab,'ylim',aa.YLim);
                            
                            goodIdx = abs(diff(b_rewardDelivered(1:idx)-t_rewardDelivered(1:idx)))<0.1;
                            b_rewardDelivered = b_rewardDelivered(goodIdx);
                            t_rewardDelivered = t_rewardDelivered(goodIdx);
                            
                            subplot(3,1,3);
                            plot(b_rewardDelivered,t_rewardDelivered,'ro');
                            hold on;
                            
                            
                            beta = [ones(length(b_rewardDelivered),1) b_rewardDelivered]\t_rewardDelivered; %Regression
                            
                            t_rewardDelivered_hat = [ones(length(b_rewardDelivered),1) b_rewardDelivered]*beta;
                            plot(b_rewardDelivered, t_rewardDelivered_hat, 'r--')
                            
                            t_stimOnset = [ones(length(b_stimOnset),1) b_stimOnset]*beta;
                            t_goCue = [ones(length(b_goCue),1) b_goCue]*beta;
                            t_responseMade = [ones(length(b_responseMade),1) b_responseMade]*beta;
                            
                            
                            %PHOTODIODE ALIGNMENT DOESNT SEEM NECESSARY
                            %                         %Then get the time of the first photodiode event
                            %                         %immediately after those t_stimOnset events, as
                            %                         %this will be the actual time the stimulus appears
                            %                         %on screen (measured by a photodiode)
                            %                         pDiode = tl.Timeline.rawDAQData(:,find(strcmp({tl.Timeline.hw.inputs.name},'photoDiode')));
                            %
                            %                         for trial = 1:length(t_stimOnset)
                            %                             %Get time idx of that stimOnset
                            %                             idx = find(t>t_stimOnset(trial),1,'first') - 1;
                            %
                            %                         end
                            
                            %CORRECT TO BEEP TIME
                            %                            audio = tl.Timeline.rawDAQData(:,find(strcmp({tl.Timeline.hw.inputs.name},'audioMonitor')));
                            %                            for trial = 1:length(t_goCue)
                            %                                %Get time idx of that stimOnset
                            %                                idx = find(t>t_goCue(trial),1,'first') - 1;
                            %                                latency_idx = find(audio(idx:end)>0.3,1,'first');
                            %                                t_goCue(trial) = t_goCue(trial) + t(latency_idx);
                            %                            end
                            
                            %                         figure;
                            %                         subplot(1,2,1); title({'RewardDelivery timebase regression',beta});
                            %                         plot(b_rewardDelivered,t_rewardDelivered,'ko',b_rewardDelivered,[ones(length(b_rewardDelivered),1) b_rewardDelivered]*beta,'k-')
                            %                         xlabel('Block timing'); ylabel('Timeline timing');
                            %                         subplot(1,2,2);
                            %                         plot(b_stimOnset,t_stimOnset,'k-')
                            %                         xlabel('Block timing'); ylabel('Timeline timing');
                            %                         title('Stimulus Onset realignment');
                            
                            behavT = [];
                            
                            behavT.timestamps = {'stimOnset','responseMade'};
                            behavT.stimOnset = t_stimOnset;
                            %                         behavT.goCue = t_goCue;
                            behavT.responseMade = t_responseMade;
                            
                            
                            %Add time events
                            
                            
                            %
                            %                         events = {tl.Timeline.hw.inputs.name};
                            %                         figure;
                            %                         for e = 1:length(events)
                            %                             subplot(5,4,e);
                            %                             pltData = tl.Timeline.rawDAQData(:, find(strcmp({tl.Timeline.hw.inputs.name},events{e} )));
                            %
                            %
                            %                             plot(t,pltData);
                            %                             title(events{e});
                            %                             xlim([0 10]);
                            %                         end
                            %
                            
                            otherInfo = { expRefs{n}{s} };
                            varargout = {W,behav,behavT,otherInfo};
                            
                            return;
                        end
                    end
                catch
                end
                
            end
            
            
        end
        
        
        
end

end