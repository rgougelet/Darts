function [temp_EEG, EEG, gui] = rm_chans_corr(temp_EEG, EEG, gui)
str = '';
chan_epoch_corrs = [];
while isempty(str)
	disp('Calculating channel pairwise correlation over time...')
	if isempty(chan_epoch_corrs)
	[n_chans, n_samps_per_epoch, n_epochs] = size(temp_EEG.data);
	if n_epochs < 2
		epoch_length_in_secs = 1/5;
		epoch_length_in_samps = temp_EEG.srate*epoch_length_in_secs;
		n_epochs = floor(temp_EEG.pnts/epoch_length_in_samps);
		n_samps_per_epoch = floor(temp_EEG.pnts/n_epochs);
		temp_data = temp_EEG.data(:,1:n_samps_per_epoch*n_epochs);
		temp_data = reshape(temp_data, n_chans, n_samps_per_epoch, n_epochs);
	else
		temp_data = temp_EEG.data;
	end
		for epoch_i = 1:n_epochs
			epoch_data = temp_data(:,:,epoch_i);
			epoch_data = diff(temp_data(:,:,epoch_i),[],2);
			epoch_corrs = abs(corr(epoch_data','type','spearman'));
			epoch_corrs(logical(eye(n_chans))) = NaN; % remove self-correlation
			chan_epoch_corrs(:,:,epoch_i) = epoch_corrs;
		end
	end
	figure(gui); clf;
	subplot(1,20,1:17);
	h = heatmap(squeeze(nanmean(chan_epoch_corrs, 2)),'GridVisible','off',...
		'XLabel', 'Correlation Over Time', 'YLabel', 'Channels', 'ColorBarVisible','off');
	h.ColorScaling = 'scaledcolumns';
	old_warning_state = warning('off', 'MATLAB:structOnObject');
	hs = struct(h);
	warning(old_warning_state);
	hs.XAxis.TickValues = [];
	subplot(1,20,18:20);
	barh(squeeze(nanmean(squeeze(nanmean(chan_epoch_corrs,3)))))
	set(gca,'YTick',[],'Ydir', 'reverse')
	
	prompt = ['\nInvestigate which channel for removal?'...
		'(Refer to plot; 0 returns to main menu.)\n'...
		'Input: '];
	str = input(prompt,'s');
	if str == '0'
		figure(gui); clf;
		break;
	elseif str2double(str)<=EEG.nbchan && str2double(str)>0
		% get channel of interest index
		coi_index = str2double(str);
		coi_data = temp_data(coi_index,:,:);
		coi_data = coi_data(:,:);
		X = [temp_EEG.chanlocs.X];
		Y = [temp_EEG.chanlocs.Y];
		Z = [temp_EEG.chanlocs.Z];
		X_dists = X(coi_index) - X;
		Y_dists = Y(coi_index) - Y;
		Z_dists = Z(coi_index) - Z;
		dists = sqrt(X_dists.^2+Y_dists.^2+Z_dists.^2);
		[~,nn_is] = mink(dists,7);
		figure(gui)
		clf;
		for nn_i = 2:length(nn_is)
			chan_i = nn_is(nn_i);
			subplot(3,2,nn_i-1);
			chan_data = temp_data(chan_i,:,:);
			chan_data = chan_data(:,:);
			scatter(chan_data,coi_data);
% 			qqplot(chan_data,coi_data);
			sgtitle(['Channel of interest ',...
				temp_EEG.chanlocs(str2double(str)).labels,...
				' against its nearest neighbors']);
			ylabel(temp_EEG.chanlocs(chan_i).labels);
		end
		prompt = ['\nRemove channel?'...
			'(y removes channel; anything else goes back.)\n'...
			'Input: '];
		str = input(prompt,'s');
		if str == 'y'
			chan_epoch_corrs = [];
			temp_EEG = pop_interp(temp_EEG, coi_index, 'spherical');
		end
	end
	str = '';
end
% trash
% 		bar(squeeze(nanmean(squeeze(nanmean(chan_epoch_corrs,2)),2)));
% 		bar(squeeze(mean(squeeze(min(chan_epoch_corrs,[],2)),2)));
% 		pw_one_minus_means = 1-squeeze(mean(chan_epoch_corrs,3)); % channels with high average correlation are suppressed
% 		pw_stds = squeeze(std(chan_epoch_corrs,0,3)); % channels with high variability in correlation are enhanced
%
% 		figure(gui); clf;
% % 		heatmap(pw_one_minus_means./pw_stds) % reciprocal coefficient of variation
% 		subplot(1,2,1);
% 		[~,min_mean_i] = min(squeeze(mean(chan_epoch_corrs,2)));
% 		plot(min_mean_i,'r.'); hold on;
% 		[~,min_med_i] = min(squeeze(median(chan_epoch_corrs,2)));
% 		plot(min_med_i,'.'); ylim([0, n_chans+1]);
% 		subplot(1,2,2);
% 		figure(gui); clf;
% 		;
%
% 		[~,min_min_i] = min(squeeze(min(chan_epoch_corrs,[],2)));
% 		plot(min_min_i,'.')
% temp_cont_data = [];
% 		for chan_i = 1:n_chans
% 			chan_data = reshape(temp_EEG.data(chan_i,:,:),1,[]);
% 			temp_cont_data(chan_i,:) = chan_data;
% 		end
% 		cont_corrs = abs(corr(temp_cont_data'));
% 		cont_corrs(logical(eye(n_chans))) = NaN;
% 		figure(gui); clf;
% 		bar(squeeze(nanmean(cont_corrs)));
% 		x = rand(20,100);
% 		y = x(:,1:50)';
% 		z = x(:,51:100)';
% 		bar(squeeze(mean((corr(y)+corr(z))/2)))
% 		bar(squeeze(mean(corr(x'))))
% 		%
% 		bar(squeeze(mean(squeeze(max(chan_epoch_corrs,[],2)),2)));
% 		min_chan_corr = squeeze(min(chan_epoch_corrs,[],3));
% 		med_chan_corr = squeeze(median(chan_epoch_corrs,3));
% 		heatmap(med_chan_corr);
% 		heatmap(1-min_chan_corr);
% 		% the problem is that the mean and median pairwise correlations for each
% 		% channel differ a lot. maybe min(min()) would work better.
% 		% 			med_chan_corrs = [];
% 		% 			for chan_i = 1:n_chans-1
% 		% 				chan_corr = med_chan_corr(:,chan_i);
% 		% 				chan_corr(chan_i) = []; % remove diagonal
% 		% 				med_chan_corrs(:,chan_i) = chan_corr;
% 		% 			end
% 		[~,min_mean_i] = min(squeeze(median(chan_epoch_corrs,2)));
% 		[~,max_mean_i] = max(squeeze(median(chan_epoch_corrs,2)));
% 		figure(gui); clf;
% 		c_min = hist(min_mean_i,1:62);
% 		c_max = hist(max_mean_i,1:62);
% 		bar([c_min;c_max]','stacked')
% 		median_pw_corrs = squeeze(median(chan_epoch_corrs,2));
% 		plot(median_pw_corrs')
% d_temp_data = diff(temp_data,1,2);
% 	corrs = abs(corr(d_temp_data','type','spearman'));
% 	for chan_i = 1:length(corrs)
% 		chan_epoch_corrs = corrs(:,chan_i);
% 		chan_epoch_corrs(chan_i) = []; % remove diagonal
% 		chans_corrs(:,chan_i) = chan_epoch_corrs;
% 	end
% 	[~,min_mean_i] = min(mean(chans_corrs));
% 	% 	non_min_is = 1:temp_EEG.nbchan;
% 	% 	non_min_is(min_i) = [];
% 	% 	nonmin_chan_data = temp_data(non_min_is,:);
% 	% 	min_chan_data = temp_data(min_i,:);
% 	figure(gui)
% 	subplot(3,1,1); boxplot(chans_corrs);
% 	subplot(3,1,2); hist(chans_corrs(min_mean_i,:))
% 	subplot(3,1,3); plot(temp_data(min_mean_i,1:(20*temp_EEG.srate)))
% 	prompt = ['\nRemove which channels?'...
% 		' (Refer to plot; 0 returns to main menu.)\n'...
% 		'Input: '];
% 	str = input(prompt,'s');
% 	if str == '0'
% 		clf;
% 		return;
% 	end
% 	if (0 < str2double(str)) && (str2double(str) < n_chans)
% 		disp('Channel removed')
% 	else
% 		disp('Enter a valid channel number.')
% 		continue;
% 	end
% 	figure(gui)
% 	corrs = abs(corrcoef(diff(temp_data')));
% 	max_corrs = maxk(corrs,2);
% 	max_corrs = max_corrs(2,:);
% 	heatmap(max_corrs)
% 	str = '';
% end