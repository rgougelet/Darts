init
pipe_dir = 'IIR BP 1 Hz Pass Edge - Lower Order - Notch Filters/';
load([data_dir,pipe_dir,'subject_struct.mat'])

%%
clearvars -except subjects

% load XLSX SNAP data
[num,txt,raw] = xlsread([data_dir,'behavioral_data_reduced.xlsx']);
r = struct;
headers = txt(1,:);
for k=1:numel(headers)
	for ri = 1:length(num)
		r(ri).(headers{k})=num(ri,k);
	end
end

% init labels for diff regression vars
chan_labs = {'Front','Back'}; % must match channels in prev script
freq_labs = {'Theta', 'Alpha'};
interval_labs = {'Delay', 'Pre_Throw'};

% initialize eeg columns to add to the behav columns
nans = num2cell(nan(length(r),1));
for chan_lab = chan_labs
		for interval_lab = interval_labs
			[r.([interval_lab{:},'_',chan_lab{:}])] = nans{:};
		end
end
%%
for subj = subjects
	for clust = subj.Clusters
		clust.Data.EOG.Filtered.Alpha.
% 		for epoch = clust.Epochs
% 		end
	end
end
%%
for subj_i = 1:length(subjects)
	for clust_i = 1:length(subjects(subj_i).Clusters)
		for epoch_i = 1:length(subjects(subj_i).Clusters(clust_i))
		epoch = subjects(subj_i).Clusters(clust_i).Epochs(epoch_i);
% 		wd = epoch.Whole.Data.EEG.Filtered.Alp
		end
	end
end
% 
% % get amps
% for subj_i = 1:length(subjs_to_include)
% 	subj_id = subjs_to_include{subj_i};
% 	load([data_dir,pipe_dir,subj_id,'_epochs']);
% 	load([data_dir,'EEG2SNAP/',subj_id,'_eeg2snap']);
% 	load([data_dir,'Latencies/',subj_id,'_eeg_512_latencies']);
% 	for chan_lab_i = 1:length(chan_labs)
% 		chan_lab = chan_labs{chan_lab_i};
% 		for freq_lab_i = 1:length(freq_labs)
% 			freq_lab = freq_labs{freq_lab_i};
% 			
% 			% init
% 			delay_trial_amps = [];
% 			prethrow_trial_amps = [];
% 			trial_eeg_throwtime = [];
% 			trial_eeg_delaytime = [];
% 			for epoch_i = 1:length(epochs)
% 				epoch = epochs{epoch_i};
% 				ep_data = epoch.post.(chan_lab);
% 				baseline_latency_i = 1;
% 				start_latency_i = 1+baseline_length_in_samples; % "latency" means sample index
% 				cue_latency_i = cue_event_latencies(epoch_i)-start_event_latencies(epoch_i);
% 				end_latency_i = (end_event_latencies(epoch_i)-offset_in_samples)-start_event_latencies(epoch_i);
% 				
% 				% timings, runs reduntantly for each chan/freq
% 				trial_eeg_delaytime(end+1) = (cue_latency_i-start_latency_i+1)/srate; % time from cue onset to dart release
% 				trial_eeg_throwtime(end+1) = (end_latency_i-cue_latency_i+1)/srate; % time from target onset to cue onset
% 				
% 				% get indices to split epoch into two parts
% 				delay_baseline_inds = 1:start_latency_i;
% 				delay_inds = start_latency_i:cue_latency_i; % target onset to throw cue onset
% 				prethrow_baseline_inds = (cue_latency_i+1-baseline_length_in_samples):cue_latency_i;
% 				prethrow_inds = cue_latency_i+1:end_latency_i;
% 
% 				bp_theta = iirsos.bp(ep_data,srate,[3 8],[2.75,8.25],.1,0);
% 				theta.amp = abs(hilbert(bp_theta)).^2;
% 				theta.delay.amp(epoch_i) = ...
% 					mean(theta.amp(delay_inds)-mean(theta.amp(delay_baseline_inds)));
% 				theta.prethrow.amp(epoch_i) = ...
% 					mean(theta.amp(prethrow_inds)-mean(theta.amp(prethrow_baseline_inds)));
% 
% 				bp_alpha = iirsos.bp(ep_data,srate,[8 12],[7.75,12.25],.1,0);
% 				alpha.amp = abs(hilbert(bp_alpha)).^2;
% 				alpha.delay.amp(epoch_i) = ...
% 					mean(alpha.amp(delay_inds)-mean(alpha.amp(delay_baseline_inds)));
% 				alpha.prethrow.amp(epoch_i) = ...
% 					mean(alpha.amp(prethrow_inds)-mean(alpha.amp(prethrow_baseline_inds)));
% 			end
% 			
% 			% assign amplitudes to matching variable col and trial row
% 			delay_amps = {r.(['Delay_',chan_lab,'_Theta'])}; % get whole columns from all subjects to add to
% 			prethrow_amps = {r.(['Pre_Throw_',chan_lab,'_Theta'])}; % get whole columns from all subjects to add to
% 			eeg_throwtime = {r.eeg_throwtime}; % get whole columns from all subjects to add to
% 			eeg_delaytime = {r.eeg_delaytime}; % get whole columns from all subjects to add to
% 			for trial_amp_i = 1:n_eeg_trials % add in subject specific trials
% 				delay_amps{eeg_to_snap_inds(trial_amp_i)} = delay_trial_amps(trial_amp_i);
% 				prethrow_amps{eeg_to_snap_inds(trial_amp_i)} = prethrow_trial_amps(trial_amp_i);
% 				eeg_throwtime{eeg_to_snap_inds(trial_amp_i)} = trial_eeg_throwtime(trial_amp_i);
% 				eeg_delaytime{eeg_to_snap_inds(trial_amp_i)} = trial_eeg_delaytime(trial_amp_i);
% 			end
% 		end
% 	end
% end