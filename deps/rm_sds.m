function [EEG, gui] = rm_sds(EEG, gui)
	str = '';
	if isempty(EEG.reject.rejmanual)
		EEG.reject.rejmanual = false(1,EEG.trials);
	end
	rej = false(size(EEG.reject.rejmanual));
	stds = mean(squeeze(std(EEG.data,0,2)),1);
	while isempty(str)
		plot_data = stds;
		plot_data(EEG.reject.rejmanual|rej) = NaN;
		plt(plot_data,gui)
		prompt = {...
			[num2str(sum(EEG.reject.rejmanual|rej)), ' epochs/chunks removed so far.\n']
			'Remove chunks/epochs with standard deviation greater\n'
			'than what? (Refer to plot; 0 returns to main menu.)\n'
			'Input: '};
		str = input([prompt{:}],'s');
		if str == '0'
			EEG.reject.rejmanual = EEG.reject.rejmanual|rej; 
			EEG = eeg_checkset(EEG);
			EEG.etc.marked = EEG.reject.rejmanual;
			clf;	return;
		end
		rej = stds>str2double(str);
		str = '';
	end
end
function plt(stat,gui)
	figure(gui)
	subplot(2,2,[1 2]); plot(stat); title('Epochs vs. Time');
	subplot(2,2,3); hist(stat); title('Epoch Std. Devs.');
	subplot(2,2,4); qqplot(stat); title('QQ');
end