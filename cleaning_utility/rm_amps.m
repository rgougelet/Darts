function [temp_EEG, EEG, gui] = rm_amps(temp_EEG, EEG, gui)
	str = '';
	clf;
	while isempty(str)
		[~, ~, n_epochs] = size(temp_EEG.data);
		if n_epochs > 1
			qs = quantile(temp_EEG.data,[0, 0.25, 0.5, 0.75, 1],2);
			mean_qs = squeeze(mean(qs,1));
			figure(gui);
			plot(mean_qs')
			prompt = ['\nRemove epochs with min/max less than/greater than what?'...
				'(Refer to plot; 0 returns to main menu.)\n'...
				'Input: '];
			str = input(prompt,'s');
			if str == '0'
				figure(gui); clf;
				break;
			end
			if str2double(str) < 0
				temp_EEG = pop_rejepoch(temp_EEG,mean_qs(1,:)<str2double(str),0);
			elseif str2double(str) > 0
				temp_EEG = pop_rejepoch(temp_EEG,mean_qs(5,:)>str2double(str),0);
			end
			qs = quantile(temp_EEG.data,[0, 0.25, 0.5, 0.75, 1],2);
			mean_qs = squeeze(mean(qs,1));
			figure(gui);
			plot(mean_qs')

			str = '';
		else
			disp(['EEG.data array has no epochs.',...
				' Consider splitting into chunks.',...
				newline]);
			break
		end
	end

%% trash