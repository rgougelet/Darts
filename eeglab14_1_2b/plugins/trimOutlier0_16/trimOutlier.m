% trimOutlier() - rejects 1.channels that are below or above the specified
%                 SD, and 2.datapoints that are above the specified
%                 threshold. Point spread width [ms] determines the range
%                 for rejection.  
%
% Usage:
%   >> EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, amplitudeThreshold, pointSpreadWidth);

% History:
% 06/27/2014 ver 1.4 by Makoto and Clement. Drastic speed up by removing for loop (thanks Christian!) Displays log (requested by Kathleen VanBenthem)
% 04/17/2014 ver 1.3 by Makoto. min(badPntsStart)=1, zero not allowed.
% 04/01/2014 ver 1.2 by Makoto. Check inputs.
% 03/26/2014 ver 1.1 by Makoto. Debug and simplify datapoint rejection
% 03/07/2014 ver 1.0 by Makoto. Former firstpassOutlierTrimmer redesigned and renamed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 09/06/2013 ver 1.5 by Makoto. No datapoint rejection reflects channel rejection.
% 08/05/2013 ver 1.4 by Makoto. Supported 3-D data (except for datapoint rejection)
% 06/27/2013 ver 1.3 by Makoto. Error message for inputting 3-D data.
% 06/13/2013 ver 1.2 by Makoto. Scalp topos added (when valid channel data exist)
% 05/22/2013 ver 1.1 by Makoto. Minor brush up.
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

function EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, amplitudeThreshold, pointSpreadWidth)

if ~(nargin==5)
    error('trimOutlier() requires 5 input arguments.')
end

% remove bad channels with SD
stdAllPnts  = std(EEG.data(:,:),0,2);
badChanIdx  = find((stdAllPnts < channelSdLowerBound) | (stdAllPnts > channelSdUpperBound));
if ~isempty(badChanIdx)
    badChanName = {EEG.chanlocs(badChanIdx).labels};
    EEG = pop_select(EEG, 'nochannel', badChanIdx);
    
    % display log
    disp(sprintf('\nThe following channels were removed:'))
    disp(badChanName)
else
    disp(sprintf('\nNo channel removed.'))
end

% return if 3-D
if length(size(EEG.data))==3
    disp('Epoched data detected: datapoint rejection will not be performed.')
    return
end

%% remove bad datapoints

% obtain the window size
windowSize = pointSpreadWidth; % millisecond
windowSizeInFrame = round(windowSize/(1000/EEG.srate)); % frame

% compute bad datapoints
absMinMaxAllChan = max([abs(min(EEG.data(:,:))); abs(max(EEG.data(:,:)))],[],1);
badPoints  = absMinMaxAllChan > amplitudeThreshold;

if any(badPoints)
    % expand badPoints
    badPointsExpanded = logical(conv(single(badPoints), ones(1,windowSizeInFrame), 'same'));
    
    % start Christian's impressive code
    rejectDataIntervals = reshape(find(diff([false badPointsExpanded false])),2,[])';
    rejectDataIntervals(:,2) = rejectDataIntervals(:,2)-1;
    
    % reject them
    EEG = pop_select(EEG, 'nopoint', [rejectDataIntervals(:,1) rejectDataIntervals(:,2)]);
    
    % display log
    badPointsInSec = length(find(badPointsExpanded))*1000/EEG.srate/1000; %#ok<*NASGU>
    disp(sprintf('\n%2.0fuV threshold with %2.0fms spreading rejected %2.1fsec data, added %1.0f boundaries.', amplitudeThreshold, windowSize, badPointsInSec, size(rejectDataIntervals,1)));
else
    disp('No datapoint rejected.');
end