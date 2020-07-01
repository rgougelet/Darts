function pipeline = get_clusters(pipeline,labels_of_each_cluster,labels_of_channels_in_each_cluster, eog_type)

tic
pipeline.nClusters = length(labels_of_each_cluster);
pipeline.ClusterLabels = labels_of_each_cluster;

subjects = pipeline.Subjects;
srate = pipeline.SampleRate;
for subj_i = 1:length(subjects)
	try
		EEG = pop_loadset('filename',subjects(subj_i).Set.name,'filepath',subjects(subj_i).Set.folder);
	catch
		eeglab;	close all;
		EEG = pop_loadset('filename',subjects(subj_i).Set.name,'filepath',subjects(subj_i).Set.folder);
	end
	
	% eeg rows must be channels
	clusters = struct;
	eog = struct;
	eeg = struct;
	for cluster_i = 1:length(labels_of_each_cluster)
		disp(['Getting Subject: ', num2str(subj_i),' Cluster: ' num2str(cluster_i)]);
		clusters(cluster_i).Index = cluster_i;
		clusters(cluster_i).SampleRate = srate;
		clusters(cluster_i).Label = labels_of_each_cluster{cluster_i};
		clusters(cluster_i).Channel.Labels = labels_of_channels_in_each_cluster{cluster_i};
		indices_for_channels_in_cluster = [];
		for chan_lab = labels_of_channels_in_each_cluster{cluster_i}
			indices_for_channels_in_cluster = [indices_for_channels_in_cluster, find(strcmpi(chan_lab,{EEG.chanlocs.labels}))];
		end
		clusters(cluster_i).Channel.Indices = indices_for_channels_in_cluster;
		
		% eog
		eog.r = get_bipolar_eog(EEG, eog_type); %% bipolar vs unipolar
		eog.f.t = iirsos.bp(eog.r,EEG.srate,[3 8],[2.75,8.25],.1,0,0);
		eog.f.a = iirsos.bp(eog.r,EEG.srate,[8 12],[7.75,12.25],.1,0,0);
		
		% eeg corrected=>filtered
		eeg.r = mean(EEG.data(indices_for_channels_in_cluster,:),1);
		eeg.c.r = eog_regression(eeg.r,eog.r);
		eeg.c.f.t = iirsos.bp(eeg.c.r,EEG.srate,[3 8],[2.75,8.25],.1,0,0);
		eeg.c.f.a = iirsos.bp(eeg.c.r,EEG.srate,[8 12],[7.75,12.25],.1,0,0);
		
		% eeg filtered=>corrected
		eeg.f.t.r = iirsos.bp(eeg.r,EEG.srate,[3 8],[2.75,8.25],.1,0,0);
		eeg.f.a.r = iirsos.bp(eeg.r,EEG.srate,[8 12],[7.75,12.25],.1,0,0);
		eeg.f.t.c = eog_regression(eeg.f.t.r, eog.f.t);
		eeg.f.a.c = eog_regression(eeg.f.a.r, eog.f.a);
		
		% map to output
		clusters(cluster_i).Data.EOG.Raw = eog.r;
		clusters(cluster_i).Data.EOG.Filtered.Theta.Raw = eog.f.t;
		clusters(cluster_i).Data.EOG.Filtered.Alpha.Raw = eog.f.a;
		
		clusters(cluster_i).Data.EEG.Raw = eeg.r;
		clusters(cluster_i).Data.EEG.Corrected.Raw = eeg.c.r;
		clusters(cluster_i).Data.EEG.Corrected.Filtered.Theta = eeg.c.f.t;
		clusters(cluster_i).Data.EEG.Corrected.Filtered.Alpha = eeg.c.f.a;
		
		clusters(cluster_i).Data.EEG.Filtered.Theta.Raw = eeg.f.t.r;
		clusters(cluster_i).Data.EEG.Filtered.Alpha.Raw = eeg.f.a.r;
		clusters(cluster_i).Data.EEG.Filtered.Theta.Corrected = eeg.f.t.c;
		clusters(cluster_i).Data.EEG.Filtered.Alpha.Corrected = eeg.f.a.c;
	end
	subjects(subj_i).Clusters = clusters;
end
pipeline.Subjects = subjects;
toc



