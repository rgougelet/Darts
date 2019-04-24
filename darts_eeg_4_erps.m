clear; close all; clc;
script_dir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis\';
cd(script_dir);
load('.\trial_data');

%% erps 
srate = 64;
epoch_half_width = 3; % value of 3 would be -3 pre-event to 3 post-event, must be symmetric

erp_t = linspace(-epoch_half_width*1000,epoch_half_width*1000,srate*epoch_half_width); % in ms
erp_baseline = [0 200]; % in ms
erp_poststim_end = 3000; % in ms

[~, erp_baseline_start_i] = min(abs(erp_t-erp_baseline(1)));
[~, erp_baseline_end_i] = min(abs(erp_t-erp_baseline(2)));
[~, erp_poststim_end_i] = min(abs(erp_t-erp_poststim_end));
erp_plot_t = erp_t(erp_baseline_start_i:erp_poststim_end_i);

target_absent_baselined = target_absent(:,erp_baseline_start_i:erp_poststim_end_i,:)-mean(target_absent(:,erp_baseline_start_i:erp_baseline_end_i,:),2);
target_absent_baselined = mean(target_absent_baselined,3)';
target_present_baselined = target_present(:,erp_baseline_start_i:erp_poststim_end_i,:)-mean(target_present(:,erp_baseline_start_i:erp_baseline_end_i,:),2);
target_present_baselined = mean(target_present_baselined,3)';
ta__tp = target_absent_baselined-target_present_baselined;

chan_labs = {'Fz','Cz','Pz','Oz'};
fig1 = figure('Position', [100 100 1024 786]);
% text(0.5, 0.98,'Target Offset; 9 subjects')
for chan_i = 1:4
	subplot(4,1,chan_i);
	plot(erp_plot_t,target_absent_baselined(:,chan_i)); hold on;
	plot(erp_plot_t,target_present_baselined(:,chan_i)); hold on;
% 	plot(plot_t,ta__tp(:,chan_i));
	legend({'Target Absent', 'Target Present'})
	title(chan_labs{chan_i})
end
saveas(fig1, 'Target Offset ERP', 'tiffn')