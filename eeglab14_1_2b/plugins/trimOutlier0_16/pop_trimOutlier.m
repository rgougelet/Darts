% pop_trimOutlier - deploys user-interactive windows showing statistical
%                   plots. Calls trimOutlier. In the line plots, black,
%                   mean across all channels; red, 2SD across all channels;
%                   blue, envelope across all channels.

% History:
% 06/27/2014 ver 1.5 by Makoto and Clement. Drastic speed up by removing for loop (thanks Christian!)
% 04/17/2014 ver 1.4 by Makoto. Channel rejection interactive part fixed. 
% 04/02/2014 ver 1.3 by Simon Due Kamronn and Makoto. Fixed com = sprintf(... %s) to %d (thanks Simon again!)
% 04/01/2014 ver 1.2 by Simon Due Kamronn and Makoto. Debugged for epoched data (thanks Simon!)
% 03/20/2014 ver 1.1 by Frank Preston and Makoto. 'windowSize = 0;' added (thanks Frank!)
% 03/07/2014 ver 1.0 by Makoto. Former firstpassOutlierTrimmer redesigned and renamed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 05/16/2013 ver 1.0 by Makoto. Created.

% Author: Makoto Miyakoshi, SCCN,INC,UCSD 2013
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG,com] = pop_trimOutlier(EEG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Does EEG.times exist? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(EEG.times)
    disp('Warning: EEG.times does not exist: Created.')
    EEG.times = (1:EEG.pnts)*(1000/EEG.srate);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Is EEG.chanlocs valid? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch ~isempty(EEG.chanlocs)
    case 0
        validChanFlag = 0;
    case 1
        switch ~isempty(EEG.chanlocs(1,1).X)
            case 0
                validChanFlag = 0;
            case 1
                validChanFlag = 1;
        end
end
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First, the topograph and bar graph for standard diviation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open main figure
handle1 = figure;
set(gcf, 'Name', 'Original Data; trimOutlier()', 'NumberTitle', 'off','Color', [0.93 0.96 1])

% compute
stdAllPnts = std(EEG.data(:,:),0,2);

if validChanFlag == 1
    % plot topograph - all channels
    subplot(3,2,1)
    topoplot(stdAllPnts, EEG.chanlocs, 'electrodes', 'on', 'emarker', {'.','k',14,1});
    colorbar
    caxis([min(stdAllPnts) max(stdAllPnts)]);
    title('Channel SD topography','FontSize', 16)

    % plot topograph - 75 percentile
    subplot(3,2,2)
    sortedPnts = sort(stdAllPnts);
    cutoffPoint = round(length(stdAllPnts)*0.75);
    cutoffValue = sortedPnts(cutoffPoint);
    std75percentIdx   = find(stdAllPnts<=cutoffValue);
    std75percentPnts  = stdAllPnts(std75percentIdx);
    chanlocs75percent = EEG.chanlocs(std75percentIdx);
    topoplot(std75percentPnts, chanlocs75percent, 'electrodes', 'on', 'emarker', {'.','k',14,1});
    colorbar
    caxis([min(std75percentPnts) max(std75percentPnts)]);
    title('75 percentile','FontSize', 16)
end

% plot bargraph
switch validChanFlag
    case 0
        subplot(2,1,1)
    case 1
        subplot(3,1,2)
end
bar(stdAllPnts);
gca;

% annotations
title('Standard deviation of all channels','FontSize', 16)
xlim([0.5 size(EEG.data,1)+0.5]);
barXLabel = get(gca,'XLabel');
barYLabel = get(gca,'YLabel');
set(barXLabel, 'String', 'Channels', 'FontSize', 16)
set(barYLabel, 'String', 'Amplitudes (uV)', 'FontSize', 16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second, time series of channels  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute
meanAllChan = mean(EEG.data(:,:));
stdAllChan  = std( EEG.data(:,:),0,1);
minAllChan  = min( EEG.data(:,:),[],1);
maxAllChan  = max( EEG.data(:,:),[],1);
posi2SDChan = meanAllChan + 2*stdAllChan;
nega2SDChan = meanAllChan - 2*stdAllChan;

% plot
switch validChanFlag
    case 0
        subplot(2,1,2)
    case 1
        subplot(3,1,3)
end
hold on;
if length(size(EEG.data))==2
    plot(EEG.times/1000, minAllChan,  'b'); plot(EEG.times/1000, maxAllChan, 'b')
    plot(EEG.times/1000, nega2SDChan, 'r'); plot(EEG.times/1000, posi2SDChan, 'r')
    plot(EEG.times/1000, meanAllChan, 'k'); 
    timeScaleEnds = [EEG.xmin EEG.xmax];
else
    tmpTimes = 1/EEG.srate*(1:length(minAllChan));
    plot(tmpTimes, minAllChan,  'b'); plot(tmpTimes, maxAllChan, 'b')
    plot(tmpTimes, nega2SDChan, 'r'); plot(tmpTimes, posi2SDChan, 'r')
    plot(tmpTimes, meanAllChan, 'k'); 
    timeScaleEnds = [0 (1/EEG.srate)*(length(minAllChan))];
end
gca;

% annotations
xlim(timeScaleEnds)
title('Time series of channel mean (Black), +/-2SD (Red), and envelope (Blue)','FontSize', 16)
lineXLabel = get(gca,'XLabel');
lineYLabel = get(gca,'YLabel');
set(lineXLabel, 'String', 'Latency (s)', 'FontSize', 16)
set(lineYLabel, 'String', 'Amplitude (uV)', 'FontSize', 16)

%%%%%%%%%%%%%%%%%%%%%%
%%% cleaning data? %%%
%%%%%%%%%%%%%%%%%%%%%%
userInput = questdlg('Proceed to clean data?','title');
if strcmp(userInput, 'No') || strcmp(userInput, 'Cancel')
    error('User canceled operation')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain channels with too large SD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
userInput = inputgui('title', 'trimOutlier()', 'geom', ...
   {{2 0 [0 0] [1 1]} {3 0 [2 0] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' 'Chanel SD upper bound'} {'style' 'edit' 'string' ''}});

if isempty(userInput{1,1})
    userInput{1,1} = Inf;
end

% open new fiugre and replicate bar graph
handle2 = figure;
set(gcf, 'Name', 'Channel rejection; trimOutlier()', 'NumberTitle', 'off', 'Color', [0.93 0.96 1])
subplot(3,1,1)
bar(stdAllPnts);
h6 = gca;

close(handle1)

% annotations
title('Before channel rejection (Red, upper bound)','FontSize', 16)
xlim([0.5 size(EEG.data(:,:),1)+0.5]);
barXLabel = get(gca,'XLabel');
barYLabel = get(gca,'YLabel');
set(barXLabel, 'String', 'Channels', 'FontSize', 16)
set(barYLabel, 'String', 'Amplitudes (uV)', 'FontSize', 16)


while ~isempty(userInput{1,1})
    if isnumeric(userInput{1,1})
        threshBar = userInput{1,1};
    else
        threshBar = str2num(userInput{1,1}); %#ok<*ST2NM>
    end
    
    % draw a threshold line
    axes(h6) %#ok<*LAXES>
    hold on
    h7 = plot(0.5:0.1:size(EEG.data(:,:),1)+0.5, threshBar, 'r');
    
    % compute
    badChans  = find(stdAllPnts > threshBar);
    goodChans = setdiff(1:EEG.nbchan, badChans);
    stdAllPntsPost = stdAllPnts(goodChans);
    disp([num2str(length(badChans)) ' channels will be rejected for upper bound threshold.'])

    % plot badChan-removed on bottom 
    subplot(3,1,2)
    bar(stdAllPntsPost);
    gca;

    % annotations
    title('After channel rejection','FontSize', 16)
    xlim([0.5 length(stdAllPntsPost)+0.5]);
    ylim([0 threshBar]);
    barXLabel = get(gca,'XLabel');
    barYLabel = get(gca,'YLabel');
    set(barXLabel, 'String', 'Channels', 'FontSize', 16)
    set(barYLabel, 'String', 'Amplitudes (uV)', 'FontSize', 16)
    
    userInput = inputgui('title', 'trimOutlier()', 'geom', ...
       {{1 1 [0 0] [1 1]}, ...
        {2 1 [0 1] [1 1]} {2 1 [1 1] [1 1]}}, ... 
    'uilist',...
       {{'style' 'text' 'string' ['Upper bound at ' num2str(threshBar) 'uV rejects ' num2str(length(badChans)) ' channels. Is it ok?']} ...
        {'style' 'text' 'string' 'If not, enter other values'} {'style' 'edit' 'string' ''}});
    if ~isempty(userInput{1,1})
        delete(h7)
    end
end
channelSdUpperBound = threshBar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain channels with too small SD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort bars
stdAllPntsPostSort = sort(stdAllPntsPost,'ascend');

% plot bars
subplot(3,1,3)
bar(stdAllPntsPostSort(1:20)); % show 20 channels from the lowest
h8 = gca;

% annotations
title('Rejecting abnormally small channels (Red, lower bound)','FontSize', 16)
xlim([0.5 20+0.5]);
barXLabel = get(gca,'XLabel');
barYLabel = get(gca,'YLabel');
set(barXLabel, 'String', 'Sorted channels', 'FontSize', 16)
set(barYLabel, 'String', 'Amplitudes (uV)', 'FontSize', 16)

userInput = inputgui('title', 'trimOutlier()', 'geom', ...
   {{2 0 [0 0] [1 1]} {3 0 [2 0] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' 'Chanel SD lower bound'} {'style' 'edit' 'string' ''}});

if isempty(userInput{1,1})
    goodChansPost = goodChans;
    threshBar = 0;
else
    axes(h8)
    while ~isempty(userInput{1,1})
        threshBar = str2num(userInput{1,1});
        
        % draw a threshold line
        hold on
        h9 = plot(0.5:0.1:size(EEG.data(:,:),1)+0.5, threshBar, 'r');
        
        % compute
        badChans      = find(stdAllPntsPost < threshBar);
        goodChansPost = setdiff(goodChans, badChans);
        disp([num2str(length(badChans)) ' channels will be rejected for lower bound threshold.'])
        
        userInput = inputgui('title', 'trimOutlier()', 'geom', ...
            {{1 1 [0 0] [1 1]}, ...
            {2 1 [0 1] [1 1]} {2 1 [1 1] [1 1]}}, ...
            'uilist',...
            {{'style' 'text' 'string' ['Lower bound at ' num2str(threshBar) 'uV rejects ' num2str(length(badChans)) ' channels. Is it ok?']} ...
            {'style' 'text' 'string' 'If not, enter other values'} {'style' 'edit' 'string' ''}});
        if ~isempty(userInput{1,1})
            delete(h9)
        end
    end
end
channelSdLowerBound = threshBar;
if channelSdLowerBound == 0
    channelSdLowerBound = -Inf;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if data not continuous, return %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(size(EEG.data))==3
    pointSpreadWidth = -Inf;
    disp('Dapapoint rejection skipped because data is not continuous.')
    EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, Inf, pointSpreadWidth);
    
    % eegh before terminate
    com = sprintf('EEG = trimOutlier(EEG, %d, %d, %d, %d);', channelSdLowerBound, channelSdUpperBound, Inf, pointSpreadWidth);
    EEG = eegh(com, EEG);
    close(handle2)
    return
end
close(handle2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reject datapoint outliers %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute again
postEEG.data    = EEG.data(goodChansPost,:);

meanAllChanPost = mean(postEEG.data);
stdAllChanPost  = std(postEEG.data,0,1);
posi2SDChanPost = meanAllChanPost + 2*stdAllChanPost;
nega2SDChanPost = meanAllChanPost - 2*stdAllChanPost;
minAllChanPost  = min(postEEG.data,[],1);
maxAllChanPost  = max(postEEG.data,[],1);

% plot time-series data
handle3 = figure;
set(gcf, 'Name', 'Datapoint Rejection; trimOutlier()', 'NumberTitle', 'off', 'Color', [0.93 0.96 1])
subplot(2,1,1)
plot(EEG.times/1000, meanAllChanPost, 'k'); hold on;
plot(EEG.times/1000, minAllChanPost,  'b'); plot(EEG.times/1000, maxAllChanPost, 'b')
plot(EEG.times/1000, nega2SDChanPost, 'r'); plot(EEG.times/1000, posi2SDChanPost, 'r')
gca;

% annotations
xlim([EEG.xmin EEG.xmax])
title('Before datapoint rejection','FontSize', 16)
lineXLabel = get(gca,'XLabel');
lineYLabel = get(gca,'YLabel');
set(lineXLabel, 'String', 'Latency (s)', 'FontSize', 16)
set(lineYLabel, 'String', 'Amplitude (uV)', 'FontSize', 16)

userInput = inputgui('title', 'trimOutlier()', 'geom', ...
   {{2 1 [0 0] [1 1]} {3 1 [2 0] [1 1]}...
    {2 1 [0 1] [1 1]} {3 1 [2 1] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' 'Enter +/- threshold [uV] (if no need, press ok)'} {'style' 'edit' 'string' ''}...
    {'style' 'text' 'string' 'Point spread width for rejection [ms]'}           {'style' 'edit' 'string' ''}});

if isempty(userInput{1,1})
    threshPnts = Inf;
    pointSpreadWidth = 0;
    windowSize = 0;
else
    while ~isempty(userInput{1,1})
        threshPnts = str2num(userInput{1,1});
        
        % obtain the window size
        windowSize = str2num(userInput{1,2}); % millisecond
        windowSizeInFrame = round(windowSize/(1000/EEG.srate)); % frame
        
        % compute bad datapoints
        absMinMaxAllChan = max([abs(minAllChanPost); abs(maxAllChanPost)],[],1);
        badPoints  = absMinMaxAllChan > threshPnts;
                
        % expand badPoints
        badPointsExpanded = logical(conv(single(badPoints), ones(1,windowSizeInFrame), 'same'));
        
        % start Christian's impressive code
        rejectDataIntervals = reshape(find(diff([false badPointsExpanded false])),2,[])';
        rejectDataIntervals(:,2) = rejectDataIntervals(:,2)-1;

        % plot how much data will be rejected in sec
        badPointsInSec = length(find(badPointsExpanded))*1000/EEG.srate/1000;
        sprintf('%2.1f sec of data will be rejected.', badPointsInSec);
        sprintf('%1.0f boundary will be made.', size(rejectDataIntervals,1));

        goodPoints = setdiff(1:EEG.pnts, find(badPointsExpanded));
        meanAllChanPostPost = meanAllChanPost(goodPoints);
        posi2SDChanPostPost = posi2SDChanPost(goodPoints);
        nega2SDChanPostPost = nega2SDChanPost(goodPoints);
        minAllChanPostPost  = minAllChanPost(goodPoints);
        maxAllChanPostPost  = maxAllChanPost(goodPoints);
        
        % plot badChan-removed on bottom
        subplot(2,1,2)
        timesPost = EEG.times(1:length(goodPoints))/1000;
        hold on;
        plot(timesPost, minAllChanPostPost,  'b'); plot(timesPost, maxAllChanPostPost, 'b')
        plot(timesPost, nega2SDChanPostPost, 'r'); plot(timesPost, posi2SDChanPostPost, 'r')
        plot(timesPost, meanAllChanPostPost, 'k');
        ylim([-threshPnts threshPnts]);
        h13 = gca;
        
        % annotations
        title('After datapoint rejection','FontSize', 16)
        xlim([0 length(timesPost)*1000/EEG.srate/1000]);
        barXLabel = get(gca,'XLabel');
        barYLabel = get(gca,'YLabel');
        set(barXLabel, 'String', 'Latency (s)', 'FontSize', 16)
        set(barYLabel, 'String', 'Amplitude (uV)', 'FontSize', 16)
        
        userInput = inputgui('title', 'trimOutlier()', 'geom', ...
            {{1 1 [0 0] [1 1]}, ...
             {3 1 [0 1] [1 1]} {3 1 [1 1] [1 1]},...
             {3 1 [0 2] [1 1]} {3 1 [1 2] [1 1]}},...
            'uilist',...
            {{'style' 'text' 'string' sprintf('Threshold %2.0fuV point spread %2.0fms rejects %2.1fsec, creates %1.0f boundaries. Is it ok?',threshPnts,windowSize,badPointsInSec,size(rejectDataIntervals,1))} ...
             {'style' 'text' 'string' 'If not, enter threshold [uV]'} {'style' 'edit' 'string' ''}...
             {'style' 'text' 'string' 'Point spread width [ms]'} {'style' 'edit' 'string' ''}});
        if ~isempty(userInput{1,1})
            cla(h13)
        end
    end
end
close(handle3)

% reject data using user-defined thresholds
EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, threshPnts, windowSize);

% eegh before terminate
com = sprintf('EEG = trimOutlier(EEG, %s, %s, %s, %s);', num2str(channelSdLowerBound), num2str(channelSdUpperBound), num2str(threshPnts), num2str(windowSize));
EEG = eegh(com, EEG);