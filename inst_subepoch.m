function subepoch_inst = inst_subepoch(cluster,start,stop)
			subepoch_inst.FirstIndex = start;
			subepoch_inst.LastIndex = stop;
			subepoch_inst.Indices = start:stop;
			subepoch_inst.Length.Samples = length(subepoch_inst.Indices); 
			subepoch_inst.Length.Seconds = subepoch_inst.Length.Samples/cluster.SampleRate;

			subepoch_inst.Data.EEG.Raw = cluster.Data.EEG.Raw(subepoch_inst.Indices);
			subepoch_inst.Data.EEG.Corrected.Raw = cluster.Data.EEG.Corrected.Raw(subepoch_inst.Indices);
			subepoch_inst.Data.EEG.Filtered.Theta.Raw = cluster.Data.EEG.Filtered.Theta.Raw(subepoch_inst.Indices);
			subepoch_inst.Data.EEG.Filtered.Alpha.Raw = cluster.Data.EEG.Filtered.Alpha.Raw(subepoch_inst.Indices);
			subepoch_inst.Data.EEG.Filtered.Theta.Corrected = cluster.Data.EEG.Filtered.Theta.Corrected(subepoch_inst.Indices);
			subepoch_inst.Data.EEG.Filtered.Alpha.Corrected = cluster.Data.EEG.Filtered.Alpha.Corrected(subepoch_inst.Indices);

			subepoch_inst.Data.EOG.Raw = cluster.Data.EOG.Raw(subepoch_inst.Indices);
			subepoch_inst.Data.EOG.Filtered.Theta.Raw = cluster.Data.EOG.Filtered.Theta.Raw(subepoch_inst.Indices);
			subepoch_inst.Data.EOG.Filtered.Alpha.Raw = cluster.Data.EOG.Filtered.Alpha.Raw(subepoch_inst.Indices);
			subepoch_inst.Data.EOG.Filtered.Theta.Corrected = cluster.Data.EOG.Filtered.Theta.Corrected(subepoch_inst.Indices);
			subepoch_inst.Data.EOG.Filtered.Alpha.Corrected = cluster.Data.EOG.Filtered.Alpha.Corrected(subepoch_inst.Indices);
end