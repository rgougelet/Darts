function [temp_EEG, EEG, gui] = split(temp_EEG, EEG, gui)
		[n_chans, ~, n_epochs] = size(temp_EEG.data);
		if n_epochs ~= 1
			disp('Your data are already epoched.')
			return
		end
		epoch_length_in_secs = 1/5;
		epoch_length_in_samps = temp_EEG.srate*epoch_length_in_secs;
		n_epochs = floor(temp_EEG.pnts/epoch_length_in_samps);
		n_samps_per_epoch = floor(temp_EEG.pnts/n_epochs);
		temp_EEG.data = temp_EEG.data(:,1:n_samps_per_epoch*n_epochs);
		temp_EEG.data = reshape(temp_EEG.data, n_chans, n_samps_per_epoch, n_epochs);
% 		new_event = temp_EEG.event;
% 		start_i = length(temp_EEG.event);
% 		for epoch_i = 1:n_epochs
% 			epoch_start_in_samps = 1+epoch_length_in_samps*epoch_i;
% 			new_event(start_i+epoch_i).type = 'cont_epoch';
% 			new_event(start_i+epoch_i).latency = epoch_start_in_samps;
% 			new_event(start_i+epoch_i).duration = 0;
% 		end
% 		temp_EEG.event = new_event;
% 		temp_EEG = pop_epoch(temp_EEG,...
% 			{'cont_epoch'},	[0 epoch_length_in_secs],...
% 			'newname', [temp_EEG.setname, '_cont_epochs'],...
% 			'epochinfo', 'yes');
		disp(['Your data are now chunked. Events not preserved.'...
			'I suggest you reject chunks only in preparation for ICA.'])
		return
end

