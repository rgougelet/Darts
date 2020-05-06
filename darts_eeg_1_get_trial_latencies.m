clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'deps/'])
data_dir = [script_dir,'data/'];
addpath(data_dir)
eeglab nogui;

%%
subjs_to_include = {
	'571'
	'579'
	'580'
	'607'
	'608'
	'616'
	'619'
	'621'
	'627'
	'631'
	};
srate = 512;

% find trial start and end times
% parfor compatible
parfor subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([data_dir,subj_id,'*',num2str(srate),'.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',data_dir);
	
	% init
	start_event_inds = [];
	start_event_latencies = [];
	start_event_strings = [];
	end_event_inds = [];
	end_event_latencies = [];
	end_event_strings = [];
	cue_event_inds = [];
	cue_event_latencies = [];
	cue_event_strings = [];
	
	% for every event marker
	for event_ind = 1:length(EEG.event)
		
		% get event type string
		event_type = num2str(EEG.event(event_ind).type);
		if length(event_type) ~= 7 % skip superfluous events
			continue
		end
		
		% find dart release event markers
		if strcmp(event_type(end-2:end),'012')
			end_event_inds(end+1) = event_ind;
			end_event_strings = [end_event_strings; event_type];
			end_event_latencies(end+1) =  round(EEG.event(event_ind).latency);
		end
	end
	
	% find prior target display onset and throw cue onset events for every dart release event
	for end_event_ind = end_event_inds
		
		% target display onset
		start_event_ind = end_event_ind;
		while start_event_ind > 0
			start_event_ind = start_event_ind-1;
			start_event_type = num2str(EEG.event(start_event_ind).type);
			if length(start_event_type) ~= 7;	continue;	end
			if strcmp(start_event_type(end-2:end),'004') ; break;	end
			if cue_event_ind==1; error('Could not find target display onset'); end
		end
		start_event_inds(end+1) = start_event_ind;
		start_event_strings = [start_event_strings; start_event_type];
		start_event_latencies(end+1) = round(EEG.event(start_event_ind).latency);
		
		% throw cue onset
		cue_event_ind = end_event_ind;
		while cue_event_ind > 0
			cue_event_ind = cue_event_ind-1;
			cue_event_type = num2str(EEG.event(cue_event_ind).type);
			if length(cue_event_type) ~= 7;	continue;	end
			if strcmp(cue_event_type(end-2:end),'009');	break; end
			if cue_event_ind==1; error('Could not find throw cue onset'); end
		end
		cue_event_inds(end+1) = cue_event_ind;
		cue_event_strings = [cue_event_strings; cue_event_type];
		cue_event_latencies(end+1) = round(EEG.event(cue_event_ind).latency);
	end
	
	% check that strings from events match for each trial
	if ~strcmp(end_event_strings(:,1:5),start_event_strings(:,1:5))
		error('Type strings dont match')
	end
	if sum(~(end_event_latencies-start_event_latencies>0))
		error('Trial starts after ends')
	end
	
	parsave([data_dir, EEG.setname,'_latencies.mat'],...
		{start_event_latencies,...
		end_event_latencies, ...
		start_event_strings, ...
		end_event_strings,...
		cue_event_strings,...
		cue_event_latencies},...
		{'start_event_latencies',...
		'end_event_latencies',...
		'start_event_strings',...
		'end_event_strings',...
		'cue_event_strings',...
		'cue_event_latencies'});
end

