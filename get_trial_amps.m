function pipeline = get_trial_amps(pipeline)
tic

%%
for subj_i = 1:pipeline.nSubjects
	subj = pipeline.Subjects(subj_i);
	for clust_i = 1:pipeline.nClusters
		clust = subj.Clusters(clust_i);
		for epoch_i = 1:subj.nEpochs
			disp(['Getting epoch amplitude: ',num2str(epoch_i)]);
			ewd = clust.Epochs(epoch_i).Whole.Data;
			cep_theta= eog_regression(...
				ewd.EEG.Filtered.Theta.Raw,...
				ewd.EOG.Filtered.Theta.Raw);
			cep_alpha = eog_regression(...
				ewd.EEG.Filtered.Alpha.Raw,...
				ewd.EOG.Filtered.Alpha.Raw);
			
			ewd.EEG.Filtered.Theta.EpochCorrected = cep_theta;
			ewd.EEG.Filtered.Alpha.EpochCorrected = cep_alpha;
			ewd.EEG.Amplitude.Theta.Raw = abs(hilbert(cep_theta));
			ewd.EEG.Amplitude.Alpha.Raw = abs(hilbert(cep_alpha));
			ewd.EEG.Amplitude.Theta.Avg = mean(ewd.EEG.Amplitude.Theta.Raw);
			ewd.EEG.Amplitude.Alpha.Avg = mean(ewd.EEG.Amplitude.Alpha.Raw);
			clust.Epochs(epoch_i).Whole.Data = ewd;
			
			ew = clust.Epochs(epoch_i).Whole;
			clust.Epochs(epoch_i).Delay.Baseline = map_whole_to_parts(ew, clust.Epochs(epoch_i).Delay.Baseline);
			clust.Epochs(epoch_i).Delay.Period = map_whole_to_parts(ew, clust.Epochs(epoch_i).Delay.Period);
			clust.Epochs(epoch_i).PreThrow.Baseline = map_whole_to_parts(ew, clust.Epochs(epoch_i).PreThrow.Baseline);
			clust.Epochs(epoch_i).PreThrow.Period = map_whole_to_parts(ew, clust.Epochs(epoch_i).PreThrow.Period);
			
			clust.Epochs(epoch_i).BaselineCorrected.Delay.Theta.Avg = ...
				clust.Epochs(epoch_i).Delay.Period.Data.EEG.Amplitude.Theta.Avg -...
				clust.Epochs(epoch_i).Delay.Baseline.Data.EEG.Amplitude.Theta.Avg;
			clust.Epochs(epoch_i).BaselineCorrected.Delay.Alpha.Avg = ...
				clust.Epochs(epoch_i).Delay.Period.Data.EEG.Amplitude.Alpha.Avg -...
				clust.Epochs(epoch_i).Delay.Baseline.Data.EEG.Amplitude.Alpha.Avg;

			clust.Epochs(epoch_i).BaselineCorrected.Delay.Theta.Raw = ...
				clust.Epochs(epoch_i).Delay.Period.Data.EEG.Amplitude.Theta.Raw -...
				clust.Epochs(epoch_i).Delay.Baseline.Data.EEG.Amplitude.Theta.Avg;
			clust.Epochs(epoch_i).BaselineCorrected.Delay.Alpha.Raw = ...
				clust.Epochs(epoch_i).Delay.Period.Data.EEG.Amplitude.Alpha.Raw -...
				clust.Epochs(epoch_i).Delay.Baseline.Data.EEG.Amplitude.Alpha.Avg;
			
			clust.Epochs(epoch_i).BaselineCorrected.PreThrow.Theta.Avg = ...
				clust.Epochs(epoch_i).PreThrow.Period.Data.EEG.Amplitude.Theta.Avg -...
				clust.Epochs(epoch_i).PreThrow.Baseline.Data.EEG.Amplitude.Theta.Avg;
			clust.Epochs(epoch_i).BaselineCorrected.PreThrow.Alpha.Avg = ...
				clust.Epochs(epoch_i).PreThrow.Period.Data.EEG.Amplitude.Alpha.Avg -...
				clust.Epochs(epoch_i).PreThrow.Baseline.Data.EEG.Amplitude.Alpha.Avg;

			clust.Epochs(epoch_i).BaselineCorrected.PreThrow.Theta.Raw = ...
				clust.Epochs(epoch_i).PreThrow.Period.Data.EEG.Amplitude.Theta.Raw -...
				clust.Epochs(epoch_i).PreThrow.Baseline.Data.EEG.Amplitude.Theta.Avg;
			clust.Epochs(epoch_i).BaselineCorrected.PreThrow.Alpha.Raw = ...
				clust.Epochs(epoch_i).PreThrow.Period.Data.EEG.Amplitude.Alpha.Raw -...
				clust.Epochs(epoch_i).PreThrow.Baseline.Data.EEG.Amplitude.Alpha.Avg;
		end
		subj.Clusters(clust_i) = clust;
	end
	pipeline.Subjects(subj.Index) = subj;
end
toc
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
% 				clust.Epochs(epoch_i) = epochs{epoch_i};
% 				ep_data = clust.Epochs(epoch_i).post.(chan_lab);
% 				baseline_latency_i = 1;
% 				start_latency_i = 1+baseline_length_in_samples; % "latency" means sample index
% 				cue_latency_i = cue_event_latencies(epoch_i)-start_event_latencies(epoch_i);
% 				end_latency_i = (end_event_latencies(epoch_i)-offset_in_samples)-start_event_latencies(epoch_i);
%

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