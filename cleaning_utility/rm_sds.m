function [temp_EEG, EEG, gui] = rm_sds(temp_EEG, EEG, gui)
str = '';
while isempty(str)
	[n_chans, n_samps_per_epoch, n_epochs] = size(temp_EEG.data);
	if n_epochs < 2
	error('Only works on epoched data')
% 		epoch_length_in_secs = 1/5;
% 		epoch_length_in_samps = temp_EEG.srate*epoch_length_in_secs;
% 		n_epochs = floor(temp_EEG.pnts/epoch_length_in_samps);
% 		n_samps_per_epoch = floor(temp_EEG.pnts/n_epochs);
% 		temp_data = temp_EEG.data(:,1:n_samps_per_epoch*n_epochs);
% 		temp_data = reshape(temp_data, n_chans, n_samps_per_epoch, n_epochs);
	end
	stds = (mean(squeeze(std(temp_EEG.data,0,2)),1));
	figure(gui)
	subplot(2,2,[1 2]); plot(stds); sgtitle('Epochs vs. Time');
	subplot(2,2,3); hist(stds); sgtitle('Epoch Std. Devs.');
	subplot(2,2,4); qqplot(stds); sgtitle('QQ');
	prompt = ...
		['\nRemove chunks with standard deviation greater',...
		' than what? (Refer to plot; 0 returns to main menu.)\n'...
		'Input: '];
	str = input(prompt,'s');
	if str == '0'
		clf
		return;
	end
	epochs_to_rej = find(stds>str2double(str));
% 	for epoch_i = 1:length(epochs_to_rej)
% 		epoch_start_i = (epochs_to_rej(epoch_i)-1)*n_samps_per_epoch+1;
% 		epoch_end_i = epochs_to_rej(epoch_i)*n_samps_per_epoch;
% 		epochs_to_rej_inds(epoch_i,:) = [epoch_start_i, epoch_end_i];
% 	end
	temp_EEG = pop_rejepoch( temp_EEG, epochs_to_rej ,0);
% 	temp_EEG = eeg_eegrej(temp_EEG, epochs_to_rej_inds);
	
% 	temp_data = temp_EEG.data;
% 	temp_data = temp_data(:,1:n_samps_per_epoch*n_epochs);
% 	temp_data = reshape(temp_data, n_chans, n_samps_per_epoch, n_epochs);
	stds = (mean(squeeze(std(temp_EEG.data,0,2)),1));
	figure(gui)
	subplot(2,2,[1 2]); plot(stds); sgtitle('Epochs vs. Time');
	subplot(2,2,3); hist(stds); sgtitle('Epoch Std. Devs.');
	subplot(2,2,4); qqplot(stds); sgtitle('QQ');
	str = '';
end
