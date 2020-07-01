function pipeline = get_epochs(pipeline,baseline_length_in_samples,offset_in_samples)
tic
subjects = pipeline.Subjects;
for subj_i = 1:length(subjects)
	lat = load([pipeline.PipeDir,'Latencies/', subjects(subj_i).ID,'_eeg_',num2str(pipeline.SampleRate),'_latencies']);
	subjects(subj_i).Index = subj_i;
	subjects(subj_i).Latencies = load([pipeline.PipeDir,'Latencies/', subjects(subj_i).ID,'_eeg_',num2str(pipeline.SampleRate),'_latencies']);
	end_latencies = lat.end_latencies;
	end_strings = lat.end_strings;
	start_latencies = lat.start_latencies;
	start_strings = lat.start_strings;
	cue_latencies = lat.cue_latencies;
	
	for cluster_i = 1:length(subjects(subj_i).Clusters)
		cluster = subjects(subj_i).Clusters(cluster_i);
		cluster.Index = cluster_i;
		epochs = [];
		% using timing data from latencies file, cut out each epoch
		for epoch_i = 1:length(start_latencies) % parfor compatible, but best not to use
			disp(['Getting epoch: ' num2str(epoch_i)]);
			epochs(epoch_i).Index = epoch_i;
			epochs(epoch_i).SubjectID = subjects(subj_i).ID;
			epochs(epoch_i).SubjectIndex = subj_i;
			epochs(epoch_i).ClusterIndex = cluster.Index;
			epochs(epoch_i).ClusterLabel = cluster.Label;
			epochs(epoch_i).DelayCondition = str2num(start_strings(epoch_i,1));
			epochs(epoch_i).TargetPosition = str2num(start_strings(epoch_i,2:3));
			epochs(epoch_i).isPres = logical(str2num(start_strings(epoch_i,4)));
			epochs(epoch_i).isAbs = ~epochs(epoch_i).isPres;
			epochs(epoch_i).TrialString = start_strings(epoch_i,1:4);
			epochs(epoch_i).SampleRate = cluster.SampleRate;
			epochs(epoch_i).Whole = inst_subepoch( cluster, ...
				start_latencies(epoch_i)-baseline_length_in_samples,...
				end_latencies(epoch_i)-offset_in_samples-1 ...
				);
			epochs(epoch_i).Delay.Baseline = inst_subepoch( cluster, ...
				start_latencies(epoch_i)-baseline_length_in_samples,...
				start_latencies(epoch_i)-1 ...
				);
			epochs(epoch_i).Delay.Period = inst_subepoch( cluster, ...
				start_latencies(epoch_i),...
				cue_latencies(epoch_i)-1 ...
				);
			epochs(epoch_i).PreThrow.Baseline = inst_subepoch( cluster, ...
				cue_latencies(epoch_i)-baseline_length_in_samples,...
				cue_latencies(epoch_i)-1 ...
				);
			epochs(epoch_i).PreThrow.Period = inst_subepoch( cluster, ...
				cue_latencies(epoch_i),...
				end_latencies(epoch_i)-offset_in_samples-1 ...
				);
		end
		subjects(subj_i).nEpochs = length(epochs);
		subjects(subj_i).Clusters(cluster_i).Epochs = epochs;
		% 			subjects(subj_i).Clusters(cluster_i).Data = [];
	end
end
pipeline.Subjects = subjects;
toc