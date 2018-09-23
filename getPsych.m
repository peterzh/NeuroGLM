function getPsych(DJSession)
% profile on

%Function which plots the psychometric curves for a session. 
% eRefs = fetchn(DJSession,'concat(session_date,"_",session_num,"_",mouse_name)->eRef');
sessions = fetch(DJSession,'*');
numSessions = length(sessions);

sessions(1).fitImg=[];
sessions(1).chronImg=[];
sessions(1).perfImg=[];
for sess = 1:numSessions
    %Get any figures matching this session
    fits = fetch(d.SessionPsychometric & sessions(sess) & 'model_name="C50-subset" OR model_name="C50-subset-2AFC"','*');
    chronos = fetch(d.SessionChronometric & sessions(sess),'*');
    perf = fetch(d.SessionPerformance & sessions(sess),'*');
    
    if ~isempty(fits)
        sessions(sess).fitImg = imread(fits.figure_path);
    end
    
    if ~isempty(chronos)
        sessions(sess).chronImg = imread(chronos.figure_path);
    end
    
    if ~isempty(perf)
        sessions(sess).perfImg = imread(perf.figure_path);
    end
end

%Sort by date


f=figure('name','',...
    'Toolbar','none','menubar','none',...
    'Position',[8 260 1904 707]);

ha_perf     = tight_subplot( 1, 5, 0.01, [0.7 0.05], 0.01);
ha_fits     = tight_subplot( 1, 5, 0.01, [0.31 0.31], 0.01);
ha_chronos  = tight_subplot( 1, 5, 0.01, [0.05 0.7], 0.01);
ha_scroll   = tight_subplot(1, 1, 0, [0 0.95], 0.02);


datenumbers = datenum({sessions.session_date});
[datenumbers,sortIdx]=sort(datenumbers);
%Sort by date
sessions = sessions(sortIdx);

% set(ha_scroll,'yticklabelmode','manual', 'TickLabelInterpreter','none');
plot(ha_scroll, datenumbers, 0, 'k.', 'markersize', 15);
xlim(ha_scroll,[min(datenumbers)-1 max(datenumbers)+1]); 


f.UserData.sessions = sessions; 
f.UserData.PlotIdx = numSessions-4;
f.UserData.subplots_perf = ha_perf;
f.UserData.subplots_fits = ha_fits;
f.UserData.subplots_chronos = ha_chronos;
f.UserData.subplot_scroll = ha_scroll;
f.UserData.datenumbers = datenumbers;


%Plot initial images
example_perfImg = sessions(find(arrayfun(@(s) ~isempty(s.perfImg), sessions),1,'first')).perfImg;
example_fitImg = sessions(find(arrayfun(@(s) ~isempty(s.fitImg), sessions),1,'first')).fitImg;
example_chronImg = sessions(find(arrayfun(@(s) ~isempty(s.chronImg), sessions),1,'first')).chronImg;
arrayfun( @(ax) imshow(example_perfImg,'parent',ax), ha_perf);
arrayfun( @(ax) imshow(example_fitImg,'parent',ax), ha_fits);
arrayfun( @(ax) imshow(example_chronImg,'parent',ax), ha_chronos);

updateDisplay(f);

set(f,'WindowScrollWheelFcn', @changeIndex);
end


function updateDisplay(figObj)

sessIdx = figObj.UserData.PlotIdx + [0:4];

%Wrap around session index
numSessions = length(figObj.UserData.sessions);
sessIdx = mod(sessIdx-1, numSessions)+1;

for i = 1:5    
    sess = figObj.UserData.sessions(sessIdx(i));
    set(get(figObj.UserData.subplots_perf(i),'children'),'CData',sess.perfImg,'UserData',sess.perfImg,'ButtonDownFcn',@zoomImg)
    set(get(figObj.UserData.subplots_fits(i),'children'),'CData',sess.fitImg,'UserData',sess.fitImg,'ButtonDownFcn',@zoomImg)
    set(get(figObj.UserData.subplots_chronos(i),'children'),'CData',sess.chronImg,'UserData',sess.chronImg,'ButtonDownFcn',@zoomImg)

    title(figObj.UserData.subplots_perf(i), sprintf('%s    %s    %d',sess.mouse_name,sess.session_date,sess.session_num) );

end

sessionDots = flipud(get(figObj.UserData.subplot_scroll, 'children'));
set(sessionDots(sessIdx),'Color',[1 0 0]);
sessionDots(sessIdx) = [];
set(sessionDots,'Color',[0 0 0]);

% drawnow;
end

function changeIndex(figObj,eventObj)
figObj.UserData.PlotIdx = figObj.UserData.PlotIdx + eventObj.VerticalScrollCount;

if figObj.UserData.PlotIdx < 1
    figObj.UserData.PlotIdx = 1;
end

if figObj.UserData.PlotIdx > ( length(figObj.UserData.sessions) - 4 )
    figObj.UserData.PlotIdx = length(figObj.UserData.sessions) - 4;
end

updateDisplay(figObj)
end

function zoomImg(obj,~)
figure('color','w','Toolbar','none','menubar','none');
imshow(obj.UserData)
end