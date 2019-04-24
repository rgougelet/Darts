% get subj 580's delaystilltimes
% initialize paths
clear; close all; clc;
scriptdir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis';
cd(scriptdir)

addpath('.\eeglab13_6_5b')
data_dir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis\data\';
addpath(data_dir)

eeglab;
close all;

EEG = pop_loadset('filename', '580_eeg_64.set','filepath',data_dir);

% for every event
start_event_inds = []; end_event_inds = [];
for event_ind = 1:length(EEG.event)
	% get event type string
	event_type = num2str(EEG.event(event_ind).type);
	if length(event_type) ~= 7
		continue
	end
	% find type strings that indicate beginning and end of delay period
	if strcmp(event_type(end-2:end),'005')
		start_event_inds = [start_event_inds EEG.event(event_ind).latency];
	end
	if strcmp(event_type(end-2:end),'010')
		end_event_inds = [end_event_inds EEG.event(event_ind).latency];
	end
end

% this variable's calculation is hardcoded into the data after
% the first subject, 580, who is also missing first 10 trials
% this variable was copy+pasted into behavioral_data.xlsx directly
delay_still_time=((end_event_inds-start_event_inds)/64)'-.1;
