rmpath('/data/common/matlab/eeglab')
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
addpath(genpath([script_dir,'deps/']));
in_dir = [script_dir, 'data/Latencies/'];
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
	} ;
srate = 512;
for subj_i = 1:length(subjs_to_include)
	subj_id = subjs_to_include{subj_i};
	%% load XLSX SNAP data
	[num,txt,raw] = xlsread([data_dir,'behavioral_data_reduced.xlsx']);
	r = struct;
	headers = txt(1,:);
	for k=1:numel(headers)
		for ri = 1:length(num)
			r(ri).(headers{k})=num(ri,k);
		end
	end
	%% correct data collection issues
	% the problem is there are some trials in the
	% xlsx file that are not in eeg
	lat = load([in_dir, subj_id,'_eeg_',num2str(srate),'_latencies']);
	end_event_latencies = lat.end_event_latencies;
	end_event_strings = lat.end_event_strings;
	start_event_latencies = lat.start_event_latencies;
	start_event_strings = lat.start_event_strings;
	cue_event_latencies = lat.cue_event_latencies;
	n_snap_trials = sum([r.subject] == str2double(subj_id));
	n_eeg_trials = length(end_event_latencies);
	eeg_trial_strs = str2num(end_event_strings(:,1:4)); % ignore warning, use str2num
	subj_inds = [r.subject] == str2double(subj_id);
	
	eeg_to_snap_inds = 1:length(eeg_trial_strs);
	if strcmp(subj_id, '580') % subj 580 missing first 10 trials
		eeg_to_snap_inds = 10+(1:n_eeg_trials);
	end
	% account for these subjects w/ missing trials
	snap_trial_strs = str2num([... % ignore warning, use str2num
		num2str([r(subj_inds).delay]'),...
		num2str([r(subj_inds).position]','%02d'),...
		num2str([r(subj_inds).pres]')]);
	if strcmp(subj_id,'616') || strcmp(subj_id,'621') || strcmp(subj_id,'627')
		eeg_to_snap_inds = [];
		for eeg_i = 1:length(eeg_trial_strs)
			for snap_i = eeg_i:length(snap_trial_strs)
				if eeg_trial_strs(eeg_i) == snap_trial_strs(snap_i)
					eeg_to_snap_inds = [eeg_to_snap_inds, snap_i];
					break
				end
			end
		end
	end
	eeg_to_snap_inds = eeg_to_snap_inds + find([r.subject]==str2double(subj_id),1) - 1;
	save([subj_id,'_eeg2snap.mat'],'eeg_to_snap_inds')
end
