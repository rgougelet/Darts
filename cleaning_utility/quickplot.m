function [temp_EEG, EEG, gui] = quickplot(temp_EEG, EEG, gui)
figure(gui);
clf;
[n_chans, n_samps_per_epoch, n_epochs] = size(temp_EEG.data);
temp_data = [];
if n_epochs > 1
		temp_data = [];
		for chan_i = 1:n_chans
			chan_data = reshape(temp_EEG.data(chan_i,:,:),1,[]);
			temp_data(chan_i,:) = chan_data;
		end
elseif n_epochs == 1
	temp_data = temp_EEG.data;
end

if ~isequaln(temp_EEG, EEG)
	figure(gui);
	subplot(2,1,1);
	qs = quantile(temp_data,[0, 0.25, 0.5, 0.75, 1],2);
	mean_qs = squeeze(mean(qs,1));
	plot(mean_qs');
	title('Cleaned data so far');
	disp('Plotted! Make sure to save your changes in the main menu.')
else
	figure(gui);
	subplot(1,1,1);
	qs = quantile(temp_EEG.data,[0, 0.25, 0.5, 0.75, 1],2);
	mean_qs = squeeze(mean(qs,1));
	plot(mean_qs');
	title('Data yet cleaned');
	disp('Plotted! You have not made any changes to the data yet.')
end


