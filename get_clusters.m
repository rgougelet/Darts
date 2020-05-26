function clusters = get_clusters(EEG, labels_of_clusters, labels_of_channels_in_each_cluster)
% eeg rows must be channels
clusters = struct;
for label_of_cluster_i = 1:length(labels_of_clusters)
	label_of_cluster = labels_of_clusters{label_of_cluster_i};
	clusters(label_of_cluster_i).Label = label_of_cluster;
	clusters(label_of_cluster_i).Channel_Labels = labels_of_channels_in_each_cluster{label_of_cluster_i};
	indices_for_channels_in_cluster = [];
	for chan_lab = labels_of_channels_in_each_cluster{label_of_cluster_i}
		indices_for_channels_in_cluster = [indices_for_channels_in_cluster, find(strcmpi(chan_lab,{EEG.chanlocs.labels}))];
	end
	
	eog.r = get_bipolar_eog(EEG);
	eog.f.t = iirsos.bp(eog.r,EEG.srate,[3 8],[2.75,8.25],.1,0);
	eog.f.a = iirsos.bp(eog.r,EEG.srate,[8 12],[7.75,12.25],.1,0);

	eeg.r = mean(EEG.data(indices_for_channels_in_cluster,:),1);
	eeg.c.r = eog_regression(eeg.r,eog.r);
	eeg.c.f.t = iirsos.bp(eeg.c.r,EEG.srate,[3 8],[2.75,8.25],.1,0);
	eeg.c.f.a = iirsos.bp(eeg.c.r,EEG.srate,[8 12],[7.75,12.25],.1,0);

	eeg.f.t.r = iirsos.bp(eeg.r,EEG.srate,[3 8],[2.75,8.25],.1,0);
	eeg.f.a.r = iirsos.bp(eeg.r,EEG.srate,[8 12],[7.75,12.25],.1,0);
	eeg.f.t.c = eog_regression(eeg.f.t.r, eog.f.t);
	eeg.f.a.c = eog_regression(eeg.f.a.r, eog.f.a);

	clusters(label_of_cluster_i).Data.EOG.Raw = eog.r;
	clusters(label_of_cluster_i).Data.EOG.Filtered.Theta = eog.f.t;
	clusters(label_of_cluster_i).Data.EOG.Filtered.Alpha = eog.f.a;

	clusters(label_of_cluster_i).Data.EEG.Raw = eeg.r;
	clusters(label_of_cluster_i).Data.EEG.Corrected.Raw = eeg.c.r;
	clusters(label_of_cluster_i).Data.EEG.Corrected.Filtered.Theta = eeg.c.f.t;
	clusters(label_of_cluster_i).Data.EEG.Corrected.Filtered.Alpha = eeg.c.f.a;
	clusters(label_of_cluster_i).Data.EEG.Filtered.Theta.Raw = eeg.f.t.r;
	clusters(label_of_cluster_i).Data.EEG.Filtered.Alpha.Raw = eeg.f.a.r;
	clusters(label_of_cluster_i).Data.EEG.Filtered.Theta.Corrected = eeg.f.t.c;
	clusters(label_of_cluster_i).Data.EEG.Filtered.Alpha.Corrected = eeg.f.a.c;

end
