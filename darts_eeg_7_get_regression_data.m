%% init
clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
addpath([script_dir,'eeglab/']);
data_dir = [script_dir,'data/'];
addpath(data_dir);
% eeglab; close;
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

% init labels for diff regression vars
chan_labs = {'Front','Back'}; % must match channels in prev script
freq_labs = {'theta', 'alpha', 'gamma'};
interval_labs = {'delay', 'pre_throw'};

%%
% load XLSX SNAP data
[num,txt,raw] = xlsread([script_dir,'behavioral_data_reduced.xlsx']);
r = struct;
headers = txt(1,:);
for k=1:numel(headers)
		for ri = 1:length(num)
			r(ri).(headers{k})=num(ri,k);
		end
end

% initialize eeg columns to add to the behav columns
nans = num2cell(nan(length(r),1));
for chan_lab = chan_labs
	for freq_lab = freq_labs
		for interval_lab = interval_labs 
			[r.([interval_lab{:},'_',chan_lab{:},'_',freq_lab{:}])] = nans{:};
		end
	end
end
[r.throwtime] = nans{:};
[r.delaytime] = nans{:};

%% not parfor compatible
for subj_i = 1:length(subjs_to_include)
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_',num2str(srate),'_clust.set'];
	
	%% correct data collection issues
	% the problem is there are some trials in the
	% xlsx file that are not in eeg
	load([data_dir, subj_id,'_eeg_',num2str(srate),'_latencies'])
	n_snap_trials = sum([r.subject] == str2double(subj_id));
	n_eeg_trials = length(end_event_latencies);
	eeg_trial_strs = str2num(end_event_strings(:,1:4)); % ignore warning, use str2num
	subj_inds = [r.subject] == str2double(subj_id);

	eeg_to_snap_inds = 1:length(eeg_trial_strs);
	if strcmp(subj_id, '580') % subj 580 missing first 10 trials
		eeg_to_snap_inds = 10+(1:n_eeg_trials);
	end
	% account for these subjects w/ missing trials
	snap_trial_strs = str2num([...
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
	
	%% filter for filter-Hilbert method
	% theta
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	EEG_theta = EEG;
	nyq = EEG_theta.srate/2;
	w0 = (5.5/nyq);
	hw = 2.5/nyq;
	wp = [w0-hw w0+hw];
	hzp = wp*nyq;
	d = .1/nyq;
	ws = [-d+wp(1) d+wp(2)];
	hzs = ws*nyq;
	rp = .01;
	rs = 6;
	[n,wn] = buttord(wp,ws,rp,rs);
	hzn = wn*nyq;
	[A,B,C,D] = butter(n,wn, 'bandpass');
	sos = ss2sos(A,B,C,D);
	x = EEG_theta.data(:,:)';
	x = sosfilt(sos,x);
	x = flip(sosfilt(sos,flip(x)));
	xx = reshape(x',size(EEG_theta.data));
	EEG_theta.data = xx;
	
	% alpha
	EEG_alpha = EEG;
	nyq = EEG_theta.srate/2;
	w0 = (10/nyq);
	hw = 2/nyq;
	wp = [w0-hw w0+hw];
	hzp = wp*nyq;
	d = .1/nyq;
	ws = [-d+wp(1) d+wp(2)];
	hzs = ws*nyq;
	rp = .01;
	rs = 6;
	[n,wn] = buttord(wp,ws,rp,rs);
	hzn = wn*nyq;
	[A,B,C,D] = butter(n,wn, 'bandpass');
	sos = ss2sos(A,B,C,D);
	x = EEG_alpha.data(:,:)';
	x = sosfilt(sos,x);
	x = flip(sosfilt(sos,flip(x)));
	xx = reshape(x',size(EEG_alpha.data));
	EEG_alpha.data = xx;
	
	% gamma
	EEG_gamma = EEG;
	nyq = EEG_gamma.srate/2;
	w0 = (142.5/nyq);
	hw = 112.5/nyq;
	wp = [w0-hw w0+hw];
	hzp = wp*nyq;
	d = .75/nyq;
	ws = [-d+wp(1) d+wp(2)];
	hzs = ws*nyq;
	rp = .01;
	rs = 6;
	[n,wn] = buttord(wp,ws,rp,rs);
	hzn = wn*nyq;
	[A,B,C,D] = butter(n,wn, 'bandpass');
	sos = ss2sos(A,B,C,D);
	x = EEG_gamma.data(:,:)';
	x = sosfilt(sos,x);
	x = flip(sosfilt(sos,flip(x)));
	xx = reshape(x',size(EEG_gamma.data));
	EEG_gamma.data = xx;
	
	% old FIR filters
	% 	EEG_theta = pop_eegfiltnew(EEG, 'locutoff', 3, 'hicutoff', 8);
	% 	EEG_alpha = pop_eegfiltnew(EEG, , 'locutoff', 8, 'hicutoff', 12);
	% 	EEG_gamma = pop_eegfiltnew(EEG, 'locutoff', 30);
	
	%% get this subjs trial-level amplitudes and timings
	EEGs = {EEG_theta, EEG_alpha, EEG_gamma};
	% at desired chans and freqs
	for chan_lab_i = 1:length(chan_labs)
		chan_lab = chan_labs{chan_lab_i};
		for freq_lab_i = 1:length(EEGs)
			freq_lab = freq_labs{freq_lab_i};
			EEG = EEGs{freq_lab_i};
			
			% hilbert transform
			EEG.data = abs(hilbert(EEG.data)).^2;
			
			% init
			delay_trial_amps = [];
			pre_throw_trial_amps = [];
			trial_throwtime = [];
			trial_delaytime = [];
			
			for eeg_trial_i = 1:n_eeg_trials
				start_latency_i = start_event_latencies(eeg_trial_i); % "latency" means sample index
				cue_latency_i = cue_event_latencies(eeg_trial_i);
				end_latency_i = end_event_latencies(eeg_trial_i)-384; % account for throw artifacts
				
				% timings, runs reduntantly for each chan/freq
				trial_throwtime(end+1) = (end_latency_i-cue_latency_i)/srate;
				trial_delaytime(end+1) = (cue_latency_i-start_latency_i)/srate;
				
				% average amplitude during delay period
				delay_baseline_inds = round(start_latency_i-0.2*srate):start_latency_i;
				delay_baseline_amp = mean(EEG.data(chan_lab_i,delay_baseline_inds),2);
				delay_inds = start_latency_i:cue_latency_i;
				delay_trial_amp = mean(EEG.data(chan_lab_i,delay_inds),2);
				delay_trial_amps(end+1) = delay_trial_amp-delay_baseline_amp;
				
				% average amplitude during pre-throw period
				pre_throw_baseline_inds = round(cue_latency_i-0.2*srate):cue_latency_i;
				pre_throw_baseline_amp = mean(EEG.data(chan_lab_i,pre_throw_baseline_inds),2);
				pre_throw_inds = cue_latency_i:end_latency_i;
				pre_throw_trial_amp = mean(EEG.data(chan_lab_i,pre_throw_inds),2);
				pre_throw_trial_amps(end+1) = pre_throw_trial_amp-pre_throw_baseline_amp;
			end
			
			% assign amplitudes to matching variable col and trial row
			delay_amps = {r.(['delay_',chan_lab,'_',freq_lab])}; % get whole columns from all subjects
			pre_throw_amps = {r.(['pre_throw_',chan_lab,'_',freq_lab])};
			throwtime = {r.throwtime};
			delaytime = {r.delaytime};
			for trial_amp_i = 1:n_eeg_trials % add in subject specific trials
				delay_amps{eeg_to_snap_inds(trial_amp_i)} = delay_trial_amps(trial_amp_i);
				pre_throw_amps{eeg_to_snap_inds(trial_amp_i)} = pre_throw_trial_amps(trial_amp_i);
				throwtime{eeg_to_snap_inds(trial_amp_i)} = trial_throwtime(trial_amp_i);
				delaytime{eeg_to_snap_inds(trial_amp_i)} = trial_delaytime(trial_amp_i);
			end
			
			% update whole columns with added subject's data
			[r.(['delay_',chan_lab,'_',freq_lab])] = delay_amps{:};
			[r.(['pre_throw_',chan_lab,'_',freq_lab])] = pre_throw_amps{:};
			[r.throwtime] = throwtime{:};
			[r.delaytime] = delaytime{:};
			
		end
	end
end
save('ic_r.mat','r')
