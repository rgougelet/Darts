init
% pipe_dir = 'IIR HP 1 Hz Pass Edge - Notch Filters/';
pipe_dir = 'IIR BP 1 Hz Pass Edge - Lower Order - Notch Filters/';
tic
for subj_i = 1:length(subjs_to_include)
	subj = [];
	subj.ID = subjs_to_include{subj_i};
	subj.Set = dir([data_dir,pipe_dir, subj.ID,'*_512.set']);
	try	
		EEG = pop_loadset('filename',subj.Set.name,'filepath',subj.Set.folder);
	catch
		eeglab;	close all;
		EEG = pop_loadset('filename',subj.Set.name,'filepath',subj.Set.folder);
	end
	
	% load trial timings
	load([data_dir,'Latencies/', subj.ID,'_eeg_',num2str(EEG.srate),'_latencies']);
	subj.Clusters = get_clusters(EEG, each_clusters_label, labels_of_channels_in_each_cluster);
	subj.SampleRate = EEG.srate;
	for cluster_i = 1:length(subj.Clusters)
		cluster = subj.Clusters(cluster_i);
		cluster.SampleRate = subj.SampleRate;

		% using timing data from latencies file, cut out each epoch
		parfor epoch_i = 1:length(start_event_latencies) % parfor compatible
			epochs(epoch_i).Whole = inst_subepoch( cluster, ...
				start_event_latencies(epoch_i)-baseline_length_in_samples,...
				end_event_latencies(epoch_i)-offset_in_samples-1 ...
				);
			epochs(epoch_i).Delay.Baseline = inst_subepoch( cluster, ...
				start_event_latencies(epoch_i)-baseline_length_in_samples,...
				start_event_latencies(epoch_i)-1 ...
				);
			epochs(epoch_i).Delay.Period = inst_subepoch( cluster, ...
				start_event_latencies(epoch_i),...
				cue_event_latencies(epoch_i)-1 ...
				);
			epochs(epoch_i).PreThrow.Baseline = inst_subepoch( cluster, ...
				cue_event_latencies(epoch_i)-baseline_length_in_samples,...
				cue_event_latencies(epoch_i)-1 ...
				);
			epochs(epoch_i).PreThrow.Period = inst_subepoch( cluster, ...
				cue_event_latencies(epoch_i),...
				end_event_latencies(epoch_i)-offset_in_samples-1 ...
				);
		end
		subj.Clusters(cluster_i).Epochs = epochs;
	end
	subjects(subj_i) = subj;
end
save([data_dir,pipe_dir,'subject_struct.mat'],'subjects', '-v7.3');
toc

%% trash
% rjgplot(x, epoch.Whole.Data.VEOG,'Seconds','uV','EOG','k'); hold on;
% 			rjgplot(x, epoch.Whole.Data.HEOG,'Seconds','uV','EOG','k');

% 			plot(epoch.Whole.Data.CorrectedRaw)
% 			plot(cluster.Data.CorrectedRaw(epoch.Whole.FirstIndex:epoch.Whole.LastIndex))
% 
% 			clf; x = 1:(epoch.Whole.Length.Samples);
% 			subplot(3,1,1);
% 			rjgplot(x, epoch.Whole.Data.CorrectedRaw,'Seconds','uV', 'Before Preprocessing');
% 			% post
% 			subplot(3,1,2);
% 			rjgplot(x,cluster.Data.CorrectedRaw(epoch.Whole.Indices),'Seconds','uV','After Preprocessing')
% 			% eog
% 			subplot(3,1,3);
% for i = 1:numel(subjects(subj_i).Clusters)

%% todo, rewrite get latency script to
%% make these names shorter.
% 	end_event_latencies = lat.end_event_latencies;
% 	end_event_strings = lat.end_event_strings;
% 	start_event_latencies = lat.start_event_latencies;
% 	start_event_strings = lat.start_event_strings;
% 	cue_event_latencies = lat.cue_event_latencies;
