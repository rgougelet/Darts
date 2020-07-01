clear; close all; clc
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath([script_dir,'eeglab/plugins/firfilt2.4/'])
addpath(genpath([script_dir,'deps/']))
data_dir = [script_dir,'data/'];
addpath(genpath(data_dir))
eeglab
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
pipe_dir = [data_dir,'IIR BP 1 Hz Pass Edge - Lower Order - Notch Filters - Laplacian/'];
%%
% find trial start and end times
% parfor compatible
parfor subj_i = 1:length(subjs_to_include)
	
	% load dataset
	subj_id = subjs_to_include{subj_i};
	subj_set = dir([pipe_dir,subj_id,'*',num2str(srate),'.set']);
	EEG = pop_loadset('filename',subj_set.name,'filepath',subj_set.folder);
	%%
	% init
	start_inds = [];
	start_latencies = [];
	start_strings = [];
	end_inds = [];
	end_latencies = [];
	end_strings = [];
	cue_inds = [];
	cue_latencies = [];
	cue_strings = [];
	
	% for every event marker
	for ind = 1:length(EEG.event)
		
		% get event type string
		type = num2str(EEG.event(ind).type);
		if length(type) ~= 7 % skip superfluous events
			continue
		end
		
		% find dart release event markers
		if strcmp(type(end-2:end),'012')
			end_inds(end+1) = ind;
			end_strings = [end_strings; type];
			end_latencies(end+1) =  round(EEG.event(ind).latency);
		end
	end
	
	% find prior target display onset and throw cue onset events for every dart release event
	for end_ind = end_inds
		
		% target display onset
		start_ind = end_ind;
		while start_ind > 0
			start_ind = start_ind-1;
			start_type = num2str(EEG.event(start_ind).type);
			if length(start_type) ~= 7;	continue;	end
			if strcmp(start_type(end-2:end),'004') ; break;	end
			if start_ind==1; error('Could not find target display onset'); end
		end
		start_inds(end+1) = start_ind;
		start_strings = [start_strings; start_type];
		start_latencies(end+1) = round(EEG.event(start_ind).latency);
		
		% throw cue onset
		cue_ind = end_ind;
		while cue_ind > start_ind
			cue_ind = cue_ind-1;
			cue_type = num2str(EEG.event(cue_ind).type);
			if length(cue_type) ~= 7;	continue;	end
			if strcmp(cue_type(end-2:end),'009');	break; end
			if cue_ind == start_ind
				error('Could not find throw cue onset before start cue');
			end
		end
		cue_inds(end+1) = cue_ind;
		cue_strings = [cue_strings; cue_type];
		cue_latencies(end+1) = round(EEG.event(cue_ind).latency);
	end
	
	% check that strings from events match for each trial
	if ~strcmp(end_strings(:,1:5),start_strings(:,1:5))
		error('Type strings dont match')
	end
	if sum(~(end_latencies-start_latencies>0))
		error('Trial starts after ends')
	end
	
	out_dir = [pipe_dir,'Latencies/'];
	mkdir(out_dir)
	parsave([out_dir, EEG.setname,'_latencies.mat'],...
		{start_latencies,...
		end_latencies, ...
		start_strings, ...
		end_strings,...
		cue_strings,...
		cue_latencies},...
		{'start_latencies',...
		'end_latencies',...
		'start_strings',...
		'end_strings',...
		'cue_strings',...
		'cue_latencies'});
end

