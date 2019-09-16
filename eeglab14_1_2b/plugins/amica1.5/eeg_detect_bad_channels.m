function [badChannels channelBadWindowRatio] = eeg_detect_bad_channels(EEG, ratioTimeBadTolerated, minAbsCorrelationThresholdAccepted)
% [badChannels channelBadWindowRatio] = eeg_detect_bad_channels(EEG, ratioTimeBadTolerated, minAbsCorrelationThresholdAccepted)
% detects bad channels and places them in badChannels output variable.
% the last two variables are optional.
% Usage:
% badChannels = eeg_detect_bad_channels(EEG);

timeWindowForCalculatingCorrelation = 1; % in seconds

if nargin<2
    ratioTimeBadTolerated = 0.01;
end;

if nargin<3
    minAbsCorrelationThresholdAccepted = 0.4;
end;

numberOfChannels = size(EEG.data,1);
numberOfFramesInTimeWindow = timeWindowForCalculatingCorrelation * EEG.srate;
numberOfTimeWindows = floor(size(EEG.data, 2) / numberOfFramesInTimeWindow);
halfNumberOfFrames = round(numberOfFramesInTimeWindow / 2);


rejectedChannels = zeros(numberOfTimeWindows, numberOfChannels);

for i=2:(numberOfTimeWindows-2) % ignore last two time windows, so we dont go out of range by chance
    eegPortion = EEG.data(:, (i*(numberOfFramesInTimeWindow-1) - halfNumberOfFrames): (i*(numberOfFramesInTimeWindow-1) + halfNumberOfFrames));
    correlationBetweenChannelsInWindow = corrcoef(eegPortion');
    correlationBetweenChannelsInWindow = correlationBetweenChannelsInWindow - diag(diag(correlationBetweenChannelsInWindow));
    absCorrelationBetweenChannelsInWindow = abs(correlationBetweenChannelsInWindow);
    maxChannelAbsCorrelationInWindow = max(absCorrelationBetweenChannelsInWindow);
    rejectedChannels(i,:)  = maxChannelAbsCorrelationInWindow < minAbsCorrelationThresholdAccepted;
end;

channelBadWindowRatio = mean(rejectedChannels,1);
badChannels = find(channelBadWindowRatio > ratioTimeBadTolerated);