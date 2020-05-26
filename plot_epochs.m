init
pipe_dir = 'IIR HP 1 Hz Pass Edge - Notch Filters/';
out_dir = [data_dir,pipe_dir,'Plot Epochs/'];
mkdir(out_dir);
overwrite = false;
for subj_i = 1:length(subjs_to_include)
	subj_id = subjs_to_include{subj_i};
	load([data_dir,pipe_dir,subj_id,'_epochs'])
	parfor epoch_i = 1:length(epochs) % parfor compatible
		epoch = epochs{epoch_i};
		% init fig
		supt = [pipe_dir(1:end-1),' ',subj_id, ' Epoch ', num2str(epoch_i)];
		
		fn = [out_dir,supt,'.png'];
		if exist(fn,'file')==2
			continue;
		else
			close all; fig = figure('Position', [0 0 1366 786], 'Visible', 'off');
			x = (1:length(epoch.inds))/srate;
			t = annotation('textbox',[.5,.98,0,0],...
				'string',supt,...
				'FitBoxToText','on',...
				'LineStyle','none',...
				'HorizontalAlignment','center',...
				'FontSize',14);
			% pre
			subplot(3,1,1);
			rjgplot(x, [epoch.pre.Front;epoch.pre.Back],'Seconds','uV');
			title('Before Preprocessing')
			legend({'Front','Back'})
			% post
			subplot(3,1,2);
			rjgplot(x, [epoch.post.Front;epoch.post.Back],'Seconds','uV');
			legend({'Front','Back'})
			title('After Preprocessing')
			% eog
			subplot(3,1,3);
			rjgplot(x, epoch.veog,'Seconds','uV','k'); hold on;
			rjgplot(x, epoch.heog,'Seconds','uV','k');
			legend({'VEOG','HEOG'})
			title('EOG')
			print(fig,fn,'-dpng');
		end
	end
end