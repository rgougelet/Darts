clear; clc; close all;
scriptdir = 'C:\Users\Rob\Desktop\darts\';
cd(scriptdir)
warning('off','MATLAB:rmpath:DirNotFound')
addpath('.\eeglab13_6_5b')
data_dir = '.\data\';
addpath(data_dir)

subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};
subjs_to_include = { '571', '608', '607', '579', '580'};

% eeglab;
close all;
% Find trial start and end times
for subj_i = 1:length(subjs_to_include)
	
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_64.set'];
	
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	
	% init
	start_event_inds = [];
	start_event_latencys = [];
	start_event_strings = [];
	end_event_inds = [];
	end_event_latencys = [];
	end_event_strings = [];
	
	% for every event
	for event_ind = 1:length(EEG.event)
		
		% get event type string
		event_type = num2str(EEG.event(event_ind).type);
		if length(event_type) ~= 7
			continue
		end
		
		% find type strings that match dart release flag 012 and store
		if strcmp(event_type(end-2:end),'012')
			end_event_inds(end+1) = event_ind;
			end_event_strings = [end_event_strings; event_type];
			end_event_latencys(end+1) =  round(EEG.event(event_ind).latency);
		end
	end
	
	% find start trial events for every dart release event
	for end_event_ind = end_event_inds
		start_event_ind = end_event_ind;
		while true
			start_event_ind = start_event_ind-1;
			start_event_type = num2str(EEG.event(start_event_ind).type);
			if length(start_event_type) ~= 7
				continue
			end
			if strcmp(start_event_type(end-2:end),'003')
				break
			end
		end
		start_event_inds(end+1) = start_event_ind;
		start_event_strings = [start_event_strings; start_event_type];
		start_event_latencys(end+1) = round(EEG.event(start_event_ind).latency);
	end
	
	% verify
	if ~strcmp(end_event_strings(:,1:5),start_event_strings(:,1:5))
		error('Type strings dont match')
	end
	if sum(~(end_event_latencys-start_event_latencys>0))
		error('Trial starts after ends')
	end
	
	parsave([data_dir, EEG.setname,'_latencys.mat'],{start_event_latencys,end_event_latencys,start_event_strings, end_event_strings},{'start_event_time_inds','end_event_time_inds','start_event_strings','end_event_strings'});
end

