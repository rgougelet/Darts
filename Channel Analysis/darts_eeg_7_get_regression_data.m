%% init
clear; close all; clc;
script_dir = '/data/common/mobi/Experiments/Darts/Analysis/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'eeglab/'])
addpath(genpath([script_dir,'deps/']))
data_dir = [script_dir,'data/'];
addpath(genpath(data_dir))

subjs_to_include = {
% 	'571'
% 	'579'
% 	'580'
% 	'607'
% 	'608'
% 	'616'
% 	'619'
	'621'
% 	'627'
% 	'631'
	};
srate = 512;
%%
cases = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 0 0; 1 0 1; 1 1 1; 1 1 0];
% not parfor compatible
% unless pwelch is fixed for parfor loops
for case_i = 1:length(cases)
	is_FIR = cases(case_i,1);
	reject_labeled_ics = cases(case_i,2);
	use_laplacian = cases(case_i,3);
	if is_FIR
		out_file = ['FIR -'];
		if reject_labeled_ics
			in_dir = ['Reject ICs -'];
			out_file = [out_file,' Reject ICs -'];
			if use_laplacian
				in_dir = [in_dir,' Laplacian'];
				out_file = [out_file,' Laplacian'];
			else
				in_dir = [in_dir,' Linked Mastoid'];
				out_file = [out_file,' Linked Mastoid'];
			end
		else
			out_file = [out_file,' Keep ICs -'];
			in_dir = ['Keep ICs -'];
			if use_laplacian
				in_dir = [in_dir,' Laplacian'];
				out_file = [out_file,' Laplacian'];
			else
				in_dir = [in_dir,' Linked Mastoid'];
				out_file = [out_file,' Linked Mastoid'];
			end
		end
	else
		out_file = ['IIR -'];
		if reject_labeled_ics
			in_dir = ['Reject ICs -'];
			out_file = [out_file,' Reject ICs -'];
			if use_laplacian
				in_dir = [in_dir,' Laplacian'];
				out_file = [out_file,' Laplacian'];
			else
				in_dir = [in_dir,' Linked Mastoid'];
				out_file = [out_file,' Linked Mastoid'];
			end
		else
			out_file = [out_file,' Keep ICs -'];
			in_dir = ['Keep ICs -'];
			if use_laplacian
				in_dir = [in_dir,' Laplacian'];
				out_file = [out_file,' Laplacian'];
			else
				in_dir = [in_dir,' Linked Mastoid'];
				out_file = [out_file,' Linked Mastoid'];
			end
		end
	end
	out_dir = [script_dir,'data/',out_file,'/']
	mkdir(out_dir);
	out_file = [out_dir,out_file];
	addpath(genpath(out_dir))
	in_dir = [script_dir,'data/',in_dir,'/']
	addpath(genpath(in_dir))

	%% init labels for diff regression vars
	chan_labs = {'Front','Back'}; % must match channels in prev script
	freq_labs = {'Theta', 'Alpha', 'Gamma'};
	interval_labs = {'Delay', 'Pre_Throw'};
	
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
	[r.eeg_throwtime] = nans{:};
	[r.eeg_delaytime] = nans{:};
	
	%% not parfor compatible
	for subj_i = 1:length(subjs_to_include)
		subj_id = subjs_to_include{subj_i};
		if use_laplacian
			ref = 'lap';
		else
			ref = 'lm';
		end
		subj_set = [subj_id,'_eeg_',num2str(srate),'_',ref,'_6_clust.set'];
		
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

		%% filter for filter-Hilbert method
		EEG = pop_loadset('filename',subj_set,'filepath',in_dir);
		
		if is_FIR
			% theta
			EEG_theta = pop_eegfiltnew(EEG, 'locutoff',3,'hicutoff',8,'plotfreqz',0);
			
			% alpha
			EEG_alpha = pop_eegfiltnew(EEG, 'locutoff',8,'hicutoff',12,'plotfreqz',0);
			
			% gamma
			EEG_gamma = pop_eegfiltnew(EEG, 'locutoff',35,'hicutoff',220,'plotfreqz',0);

		else
			% theta
			EEG_theta = EEG;
			[EEG_theta.data]= iirsos.bp(EEG_theta.data, EEG_theta.srate,[3 8], [2.75, 8.25],10e-3,0);
			
			% alpha
			EEG_alpha = EEG;
			EEG_alpha.data = iirsos.bp(EEG_alpha.data, EEG_alpha.srate,[8 12], [7.75, 12.25],10e-3,0);
			
			% gamma
			EEG_gamma = EEG;
			EEG_gamma.data = iirsos.bp(EEG_gamma.data, EEG_gamma.srate,[35 220], [30 225],10e-3,0);
		end

		%% get this subjs trial-level amplitudes and timings
		plot_EEGs = {EEG, EEG_theta, EEG_alpha, EEG_gamma};
		for plot_i = 1:length(plot_EEGs)
			EEG = plot_EEGs{plot_i};
			
			% trim and plot to verify filters
			plot_EEG = EEG;
			plot_EEG.data = [];
			for event_i = 1:length(start_event_latencies)
				epoch = EEG.data(:,start_event_latencies(event_i):end_event_latencies(event_i)-384); % 384 to correct for motion artifacts
				plot_EEG.data = [plot_EEG.data,epoch-mean(epoch,2)];
				if plot_i == 1
					figure('Visible', 'off');
					plot(epoch'-mean(epoch,2))
					t = [subj_id, ' Epoch ', num2str(event_i)];
					title(t);
					saveas(gcf,[out_dir,t,'.jpg']);
				end
			end
			figure('Visible', 'off'); pwelch(plot_EEG.data(:,:)',5000,20,[],512,'onesided');
			ylim([-30 100]);
			title(subj_id);
			saveas(gcf,[out_file,' - Post-ICA_Trimmed_', subj_id, '_', num2str(plot_i),'.jpg']);
			close
		end
		
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
				trial_eeg_throwtime = [];
				trial_eeg_delaytime = [];
				offset_in_samples = 384; % cut off end of trial to account for throw artifacts
				for eeg_trial_i = 1:n_eeg_trials
					start_latency_i = start_event_latencies(eeg_trial_i); % "latency" means sample index
					cue_latency_i = cue_event_latencies(eeg_trial_i);
					end_latency_i = end_event_latencies(eeg_trial_i);
					
					% timings, runs reduntantly for each chan/freq
					trial_eeg_throwtime(end+1) = (end_latency_i-cue_latency_i)/srate; % time from target onset to cue onset
					trial_eeg_delaytime(end+1) = (cue_latency_i-start_latency_i)/srate; % time from cue onset to dart release
					
					% average amplitude during delay period
					delay_baseline_inds = round(start_latency_i-0.2*srate):start_latency_i; % 200 ms baseline
					delay_baseline_amp = mean(EEG.data(chan_lab_i,delay_baseline_inds),2);
					delay_inds = start_latency_i:cue_latency_i; % target onset to throw cue onset
					delay_trial_amp = mean(EEG.data(chan_lab_i,delay_inds),2);
					delay_trial_amps(end+1) = delay_trial_amp-delay_baseline_amp;
					
					% average amplitude during pre-throw period
					pre_throw_baseline_inds = round(cue_latency_i-0.2*srate):cue_latency_i; % 200 ms baseline
					pre_throw_baseline_amp = mean(EEG.data(chan_lab_i,pre_throw_baseline_inds),2);
					pre_throw_inds = cue_latency_i:end_latency_i-offset_in_samples; % throw cue onset to dart release minus offset
					pre_throw_trial_amp = mean(EEG.data(chan_lab_i,pre_throw_inds),2);
					pre_throw_trial_amps(end+1) = pre_throw_trial_amp-pre_throw_baseline_amp;
				end
				
				% assign amplitudes to matching variable col and trial row
				delay_amps = {r.(['Delay_',chan_lab,'_',freq_lab])}; % get whole columns from all subjects
				pre_throw_amps = {r.(['Pre_Throw_',chan_lab,'_',freq_lab])};
				eeg_throwtime = {r.eeg_throwtime};
				eeg_delaytime = {r.eeg_delaytime};
				for trial_amp_i = 1:n_eeg_trials % add in subject specific trials
					delay_amps{eeg_to_snap_inds(trial_amp_i)} = delay_trial_amps(trial_amp_i);
					pre_throw_amps{eeg_to_snap_inds(trial_amp_i)} = pre_throw_trial_amps(trial_amp_i);
					eeg_throwtime{eeg_to_snap_inds(trial_amp_i)} = trial_eeg_throwtime(trial_amp_i);
					eeg_delaytime{eeg_to_snap_inds(trial_amp_i)} = trial_eeg_delaytime(trial_amp_i);
				end
				
				% update whole columns with added subject's data
				[r.(['Delay_',chan_lab,'_',freq_lab])] = delay_amps{:};
				[r.(['Pre_Throw_',chan_lab,'_',freq_lab])] = pre_throw_amps{:};
				[r.eeg_throwtime] = eeg_throwtime{:};
				[r.eeg_delaytime] = eeg_delaytime{:};
				
			end
		end
	end
	parsave([out_file,'.mat'],r,'r')
end
