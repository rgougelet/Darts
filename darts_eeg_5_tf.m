clear; close all; clc;
script_dir = 'C:\Users\Rob\Desktop\Dropbox\darts_eeg_analysis\';
cd(script_dir);
load('.\trial_data');

% data parameters
srate = 64;
epoch_half_width = 3;
chan_labs = {'Fz','Cz','Pz','Oz'};

% wavelet parameters
num_frex = 128;
min_freq =  1;
max_freq = 32;

% set range for variable number of wavelet cycles
range_cycles = [ 3 8 ];

% other wavelet parameters
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
wavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);
tf_time = -epoch_half_width:(1/srate):epoch_half_width;
half_wave_size = (length(tf_time)-1)/2;

% tf parameters
tf_t = linspace(-epoch_half_width*1000,epoch_half_width*1000,srate*epoch_half_width); % in ms
tf_baseline = [-100 0]; % in ms
tf_poststim_end = 250; % in ms
[~, tf_baseline_start_i] = min(abs(tf_t-tf_baseline(1)));
[~, tf_baseline_end_i] = min(abs(tf_t-tf_baseline(2)));
[~, tf_poststim_end_i] = min(abs(tf_t-tf_poststim_end));
tf_plot_t = tf_t(tf_baseline_start_i:tf_poststim_end_i);

%% morelet wavelet transform
data_labs = {'ta','tp'};

for data_i = 1:length(data_labs)
	
	datas = {target_absent,target_present};
	data_label = data_labs{data_i};
	data = datas{data_i};
	
	% FFT parameters
	nWave = length(tf_time);
	n_chans = size(data,1);
	n_trials = size(data,3);
	n_pnts = size(data,2);
	n_data = n_chans*n_trials*n_pnts;
	n_conv = nWave+n_data-1;
	
	% initialize output time-frequency data
	tf = zeros(num_frex,n_chans, n_pnts, n_trials);
	
	% now compute the FFT of all trials concatenated
	all_data = reshape(data ,1,[]);
	dataX   = fft(all_data, n_conv);
	
	% loop over frequencies
	tic
	for fi=1:num_frex
		% create wavelet and get its FFT
		s = wavecycles(fi)/(2*pi*frex(fi));
		wavelet  = exp(2*1i*pi*frex(fi).*tf_time) .* exp(-tf_time.^2./(2*s^2));
		waveletX = fft(wavelet,n_conv);
		waveletX = waveletX./max(waveletX);
		
		% now run convolution in one step
		as = ifft(waveletX .* dataX);
		as = as(half_wave_size+1:(end-half_wave_size));
		
		% and reshape back to time X trials
		as = reshape( as, n_chans, n_pnts, n_trials );
		
		% compute power and average over trials
		tf(fi,:,:,:) = abs(as);
	end
	toc
	tf_baselined = tf(:,:,tf_baseline_start_i:tf_poststim_end_i,:)-mean(tf(:,:,tf_baseline_start_i:tf_baseline_end_i,:),3);
	parsave([data_label,'_tf'],tf_baselined, [data_label,'_tf'])
end

%%
clear;
load('ta_tf.mat')
load('tp_tf.mat')

all_trials = cat(4,ta_tf,tp_tf);
n_ta = size(ta_tf,4);

tic
for perm_i = 1:1000
	rnd_trials = all_trials(:,:,:,randperm(size(all_trials,4)));
	rnd_tas(:,:,:,perm_i) = mean(rnd_trials(:,:,:,1:n_ta),4);
	rnd_tps(:,:,:,perm_i) = mean(rnd_trials(:,:,:,n_ta+1:end),4);
end
	parsave('rnd_tf',{rnd_tas,rnd_tps}, {'rnd_tas','rnd_tps'})
toc

%%
load('ta_tf.mat')
load('tp_tf.mat')
load('rnd_tf.mat')
n_ta = size(ta_tf,4);
n_tp = size(tp_tf,4);
all_trials = cat(4,ta_tf,tp_tf);
orig_ta=mean(ta_tf,4);
orig_tp=mean(tp_tf,4);
orig_diff = orig_ta-orig_tp;
% rnd_diffs = rnd_tas-rnd_tps;
figure; contourf(tf_plot_t,frex,orig_diff,100,'linecolor','none');colorbar;
figure; contourf(tf_plot_t,frex,abs(orig_diff),100,'linecolor','none');colorbar;
a = 0.05;
n_cmprsns = size(all_trials,1)*size(all_trials,2);
a_bonf = a/n_cmprsns;
extr_cnt = sum(abs(rnd_diffs(:,:,:))>abs(orig_diff),3)%/1000>(50/1000*n_cmprsns);
figure;contourf(tf_plot_t,frex,extr_cnt,100,'linecolor','none'); colorbar;


%% Plot
% Fz
load('Fz_ta_tf.mat')
load('Fz_tp_tf.mat')

fig = figure('Position', [100 100 1024 786]);
text(0.5, 0.98,'Target Onset; 10 subjects')

subplot(4,3,1);
contourf(tf_plot_t,frex,Fz_ta_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Fz Target Absent')

subplot(4,3,2);
contourf(tf_plot_t,frex,Fz_tp_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Fz Target Present')

subplot(4,3,3);
contourf(tf_plot_t,frex,Fz_ta_tf-Fz_tp_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Fz TA>TP')

% Cz
load('Cz_ta_tf.mat')
load('Cz_tp_tf.mat')

subplot(4,3,4);
contourf(tf_plot_t,frex,Cz_ta_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Cz Target Absent')

subplot(4,3,5);
contourf(tf_plot_t,frex,Cz_tp_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Cz Target Present')

subplot(4,3,6);
contourf(tf_plot_t,frex,Cz_ta_tf-Cz_tp_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Cz TA>TP')

% Pz
load('Pz_ta_tf.mat')
load('Pz_tp_tf.mat')

subplot(4,3,7);
contourf(tf_plot_t,frex,Pz_ta_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Pz Target Absent')

subplot(4,3,8);
contourf(tf_plot_t,frex,Pz_tp_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Pz Target Present')

subplot(4,3,9);
contourf(tf_plot_t,frex,Pz_ta_tf-Pz_tp_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Pz TA>TP')

% Oz
load('Oz_ta_tf.mat')
load('Oz_tp_tf.mat')

subplot(4,3,10);
contourf(tf_plot_t,frex,Oz_ta_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Oz Target Absent')

subplot(4,3,11);
contourf(tf_plot_t,frex,Oz_tp_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Oz Target Present')

subplot(4,3,12);
contourf(tf_plot_t,frex,Oz_ta_tf-Oz_tp_tf,100,'linecolor','none')
colorbar;
set(gca,'YScale', 'log')
set(gca,'ytick', ceil(logspace(log10(1),log10(max(frex)),10)))
title('Oz TA>TP')

saveas(fig, 'Target Onset TF', 'tiffn')