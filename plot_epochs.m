function plot_epochs(pipeline,overwrite)
pipe_dir = pipeline.PipeDir;
out_dir = [pipe_dir,'Plot Epochs/'];
mkdir(out_dir);
pipe_name = pipeline.Name;
nclusters = pipeline.nClusters;
for subj_i = 1:pipeline.nSubjects
	subj = pipeline.Subjects(subj_i);
	id = subj.ID;
	subj.Clusters = rmfield(subj.Clusters,'Data');
	parfor epoch_i = 1:subj.nEpochs % parfor compatible, could be optimized more tho
		disp(['Plotting epoch: ',num2str(epoch_i)]);
		% init fig
		close all;
		fig = figure('Position', [0 0 1366 786], 'Visible', 'off');
		super_title = [pipe_name,' - Subject ',id, ' Epoch ', num2str(epoch_i)];
		fn = [out_dir,super_title,'.png'];
		if exist(fn,'file')==2 && ~overwrite
			continue;
		end
		% plots both clusters
		for clust_i = 1:nclusters
			clust = subj.Clusters(clust_i);
			epoch = clust.Epochs(epoch_i);
			
			x = (1:length(epoch.Whole.Indices))/clust.SampleRate;
			annotation('textbox',[.5,.98,0,0],...
				'string',super_title,...
				'FitBoxToText','on',...
				'LineStyle','none',...
				'HorizontalAlignment','center',...
				'FontSize',14);
			ed = epoch.Whole.Data;
			ceeg = eog_regression(ed.EEG.Raw,ed.EOG.Raw);
			
			set(0, 'CurrentFigure', fig)
			subplot(4,1,1); hold on;
			plot(x, ed.EEG.Raw,'DisplayName',clust.Label); xlabel('Seconds'); ylabel('uV');
			legend show;
			title('Before Correction');
			subplot(4,1,2); hold on;
			plot(x, ed.EEG.Corrected.Raw,'DisplayName',clust.Label); xlabel('Seconds'); ylabel('uV');
			legend show;
			title('After Correction using All Data');
			subplot(4,1,3); hold on;
			plot(x, ceeg,'DisplayName',clust.Label); xlabel('Seconds'); ylabel('uV');
			legend show;
			title('After Correction using Only Epoch Data');

% 			figure; clf; 
% 			cep_theta = eog_regression(...
% 				ed.EEG.Filtered.Theta.Raw,...
% 				ed.EOG.Filtered.Theta.Raw);
% 			cep_alpha = eog_regression(...
% 				ed.EEG.Filtered.Alpha.Raw,...
% 				ed.EOG.Filtered.Alpha.Raw);
% 			subplot(2,1,1); hold on;
% 			plot(ed.EEG.Raw); 
% 			plot(ed.EEG.Filtered.Theta.Raw)
% 			plot(ed.EEG.Filtered.Alpha.Raw)
% 			title('Without Correction')
% 			legend({'Raw','Theta','Alpha'});
% 			subplot(2,1,2); hold on;
% 			plot(ceeg); 
% 			plot(cep_theta);
% 			plot(cep_alpha);
% 			title('With Correction')
% 			legend({'Raw','Theta','Alpha'});
		end
		subplot(4,1,4); hold on;
		h = plot(x, ed.EOG.Raw'); xlabel('Seconds'); ylabel('uV');
		if size(ed.EOG.Raw,1) == 2
			set(h, {'color'}, {[0 0 0]; [1 0 0]});
			legend({'VEOG','HEOG'});
		else
			set(h, {'color'}, {[0 0 0]; [1 0 0]; [1 0 1]});
			legend({'VEOG','LHEOG','RHEOG'});
		end
		print(fig,fn,'-dpng');
	end
end

