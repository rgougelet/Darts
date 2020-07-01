%% plot epoch time series
et = []; ea = [];
for subj = pipeline.Subjects
	for clust = subj.Clusters
		for epoch = clust.Epochs
			et(end+1) = length(epoch.BaselineCorrected.Delay.Theta.Raw);
			ea(end+1) = length(epoch.BaselineCorrected.Delay.Alpha.Raw);
		end
	end
end