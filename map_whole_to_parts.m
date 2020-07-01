	function e_to = map_whole_to_parts(e_from,e_to)
	inds = e_to.Indices-e_from.FirstIndex+1;
	e_to.Data.EEG.Filtered.Theta.EpochCorrected = e_from.Data.EEG.Filtered.Theta.EpochCorrected(inds);
	e_to.Data.EEG.Filtered.Alpha.EpochCorrected = e_from.Data.EEG.Filtered.Alpha.EpochCorrected(inds);
	e_to.Data.EEG.Amplitude.Theta.Raw = e_from.Data.EEG.Amplitude.Theta.Raw(inds);
	e_to.Data.EEG.Amplitude.Alpha.Raw = e_from.Data.EEG.Amplitude.Alpha.Raw(inds);
	e_to.Data.EEG.Amplitude.Theta.Avg = mean(e_to.Data.EEG.Amplitude.Theta.Raw);
	e_to.Data.EEG.Amplitude.Alpha.Avg = mean(e_to.Data.EEG.Amplitude.Alpha.Raw);
	end