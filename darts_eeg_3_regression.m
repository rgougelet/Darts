clear; close all; clc;
script_dir = '/data/mobi/Darts/Analysis/Analysis_Sept-2019/darts/';
cd(script_dir);
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/data/common/matlab/eeglab')
addpath([script_dir,'/eeglab2019_0/'])
eeglab;
close all;
data_dir = [script_dir,'/data/'];
addpath(data_dir)

subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};
srate = 512;

% init labels for diff regression vars
chan_labs = {'Fz','Cz','Pz','Oz'}; % must match channels in prev script
freq_labs = {'theta', 'alpha', 'gamma'};
interval_labs = {'delay', 'pre_throw'};

% load XLSX SNAP data
[num,txt,raw] = xlsread([script_dir,'behavioral_data_reduced.xlsx']);
headers = txt(1,:);
for k=1:numel(headers)
	xlsx.(headers{k})=num(:,k) ;
end

%%
% initialize eeg columns to add to the behav columns
for chan_lab = chan_labs
	for freq_lab = freq_labs
		for interval_lab = interval_labs
			xlsx.([interval_lab{:},'_',chan_lab{:},'_',freq_lab{:}]) = nan(length(xlsx.pres),1);
		end
	end
end
xlsx.eeg_throwtime = nan(length(xlsx.pres),1);
xlsx.eeg_delaytime = nan(length(xlsx.pres),1);

% for each subject
for subj_i = 1:length(subjs_to_include)
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_',num2str(srate),'_laplace.set'];
	
	% bandpass filter at desired freqs
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	EEG_theta = pop_eegfiltnew(EEG, 3, 7);
	EEG_alpha = pop_eegfiltnew(EEG, 8, 12);
	EEG_gamma = pop_eegfiltnew(EEG, 30);
	EEGs = {EEG_theta,EEG_alpha, EEG_gamma};
	
	% correct data collection issues
	% the problem is there are some trials in the
	% xlsx file that are not in eeg
	load([data_dir, subj_id,'_eeg_',num2str(srate),'_latencys'])
	n_snap_trials = sum(xlsx.subject == str2double(subj_id));
	n_eeg_trials = length(end_event_latencys);
	eeg_trial_strs = str2num(end_event_strings(:,1:4)); % ignore warning, use str2num
	subj_inds = xlsx.subject == str2double(subj_id);
	snap_trial_strs = str2num([num2str(xlsx.delay(subj_inds)),...
		num2str(xlsx.position(subj_inds),'%02d'),num2str(xlsx.pres(subj_inds))]);
	eeg_to_snap_inds = 1:length(eeg_trial_strs);
	if strcmp(subj_id, '580') % subj 580 missing first 10 trials
		eeg_to_snap_inds = 10+(1:n_eeg_trials);
	end
	% account for these subjects w/ missing trials
	if strcmp(subj_id,'616') || strcmp(subj_id,'621') || strcmp(subj_id,'627')
		eeg_to_snap_inds = [];
		for eeg_i = 1:length(eeg_trial_strs)
			for snap_i = eeg_i:length(snap_trial_strs)
				if eeg_trial_strs(eeg_i) == snap_trial_strs(snap_i)
					eeg_to_snap_inds = [eeg_to_snap_inds, snap_i];
					break
				end
			end
		end
	end
	eeg_to_snap_inds = eeg_to_snap_inds + find(xlsx.subject==str2double(subj_id),1) - 1;
	
	% retrieve this subjs trial-level amplitudes and timings
	% at desired chans and freqs
	for chan_lab_i = 1:length(chan_labs)
		chan_lab = chan_labs{chan_lab_i};
		for freq_lab_i = 1:length(EEGs)
			freq_lab = freq_labs{freq_lab_i};
			EEG = EEGs{freq_lab_i};
			
			% hilbert transform
			EEG.data = abs(hilbert(EEG.data)).^2;
			
			% init
			delay_trial_amps = [];
			pre_throw_trial_amps = [];
			trial_throwtime = [];
			trial_delaytime = [];
			
			for eeg_trial_i = 1:n_eeg_trials
				start_latency_i = start_event_latencys(eeg_trial_i); % "latency" means sample index
				cue_latency_i = cue_event_latencys(eeg_trial_i);
				end_latency_i = end_event_latencys(eeg_trial_i);
				
				% average amplitude during delay period
				delay_baseline_inds = round(start_latency_i-0.2*srate):start_latency_i;
				delay_baseline_amp = mean(EEG.data(chan_lab_i,delay_baseline_inds),2);
				delay_inds = start_latency_i:cue_latency_i;
				delay_trial_amp = mean(EEG.data(chan_lab_i,delay_inds),2);
				delay_trial_amps(end+1) = delay_trial_amp-delay_baseline_amp;
				
				% average amplitude during pre-throw period
				pre_throw_baseline_inds = round(cue_latency_i-0.2*srate):cue_latency_i;
				pre_throw_baseline_amp = mean(EEG.data(chan_lab_i,pre_throw_baseline_inds),2);
				pre_throw_inds = cue_latency_i:end_latency_i;
				pre_throw_trial_amp = mean(EEG.data(chan_lab_i,pre_throw_inds),2);
				pre_throw_trial_amps(end+1) = pre_throw_trial_amp-pre_throw_baseline_amp;
				
				% timings, runs reduntantly for each chan/freq
				trial_throwtime(end+1) = (end_latency_i-cue_latency_i)/srate;
				trial_delaytime(end+1) = (cue_latency_i-start_latency_i)/srate;
			end
			
			% assign amplitudes to matching variable col and trial row
			delay_amps = xlsx.(['delay_',chan_lab,'_',freq_lab]);
			pre_throw_amps = xlsx.(['pre_throw_',chan_lab,'_',freq_lab]);
			throwtime = xlsx.eeg_throwtime;
			delaytime = xlsx.eeg_delaytime;
			for trial_amp_i = 1:n_eeg_trials
				delay_amps(eeg_to_snap_inds(trial_amp_i)) = delay_trial_amps(trial_amp_i);
				pre_throw_amps(eeg_to_snap_inds(trial_amp_i)) = pre_throw_trial_amps(trial_amp_i);
				throwtime(eeg_to_snap_inds(trial_amp_i)) = trial_throwtime(trial_amp_i);
				delaytime(eeg_to_snap_inds(trial_amp_i)) = trial_delaytime(trial_amp_i);
			end
			xlsx.(['delay_',chan_lab,'_',freq_lab]) = delay_amps;
			xlsx.(['pre_throw_',chan_lab,'_',freq_lab]) = pre_throw_amps;
			xlsx.eeg_throwtime = throwtime;
			xlsx.eeg_delaytime = delaytime;
			
		end
	end
end
parsave('xlsx.mat',xlsx,'xlsx')

%% make regression vars
load('xlsx.mat');
keep_inds = ~isnan(xlsx.distance);
orig_throwtime = xlsx.eeg_throwtime(keep_inds);
orig_delaytime = xlsx.eeg_delaytime(keep_inds);
orig_distance = xlsx.distance(keep_inds);
t_distance = orig_distance.^(1/2); % sqrt(paper units)
pres = categorical(xlsx.pres(keep_inds));

% convert to z scores
delaytime = (orig_delaytime-nanmean(orig_delaytime))/nanstd(orig_delaytime);
distance = (t_distance-nanmean(t_distance))/nanstd(t_distance);
throwtime = (orig_throwtime-nanmean(orig_throwtime))/nanstd(orig_throwtime);
varnames = {'distance'; 'throwtime'; 'delaytime';'pres'};

% create workspace variables for regression
chan_labs = {'Fz','Cz','Pz','Oz'};
freq_labs = {'theta', 'alpha', 'gamma'};
interval_labs = {'delay', 'pre_throw'};
for chan_lab = chan_labs
	for freq_lab = freq_labs
		for interval_lab = interval_labs
			varname = [interval_lab{:},'_',chan_lab{:},'_',freq_lab{:}];
			varnames{end+1} = varname;
			eval([varname,' = (xlsx.',varname,'(keep_inds)-nanmean(xlsx.',varname,'(keep_inds)))/nanstd(xlsx.',varname,'(keep_inds));'])
		end
	end
end
eval( ['reg_table = table(', strjoin(varnames',','),');'] );

%% robust regression
mdl_robust = fitlm(reg_table,'interactions','ResponseVar', 'distance', 'RobustOpts','on');
summary_results = sortrows(anova(mdl_robust),5,'ascend');
coeff_results = sortrows(mdl_robust.Coefficients,4,'ascend');
save('results.mat');

% save results to xls file for table and fig 5 and 6
writetable(coeff_results,'full_model_coeffs.xls','WriteRowNames',true) 

%% results
load('results.mat');

% behavioral results
orig_mean = nanmean(orig_distance); % in original paper units
orig_sd = nanstd(orig_distance);
irl_mean = nanmean(orig_distance)*6.55 % actual distance from target to dart "in real life"
irl_median = nanmedian(orig_distance)*6.55 % actual distance from target to dart "in real life"
irl_sd = nanstd(orig_distance)*6.55

orig_throwtime_mean = nanmean(orig_throwtime)
orig_throwtime_median = nanmedian(orig_throwtime)
orig_throwtime_sd = nanstd(orig_throwtime)

sum(xlsx.pres(keep_inds))
sum(~xlsx.pres(keep_inds))

% regression results
mdl_robust
sum(mdl_robust.Robust.Weights==0)
sortrows(mdl_robust.Residuals.Raw(mdl_robust.Robust.Weights==0))

z_errors = mdl_robust.Residuals.Raw; % in z units
t_unz_errors = (z_errors+nanmean(t_distance)).*nanstd(t_distance); % transformed paper units
unz_errors = t_unz_errors.^2; % in paper units
irl_unz_errors = unz_errors*6.55; % at projection screen scale
irl_mad_unz_error = nanmedian(abs(irl_unz_errors))

% figure 2a, exclude outliers for plotting sake
figure;
hist(mdl_robust.Residuals.Raw(mdl_robust.Robust.Weights~=0),100);
xlim([-5,5]);
saveas(gcf,'fig2a.jpg')

% figure 2b
figure;
y_hat = predict(mdl_robust, reg_table);
plot(reg_table.distance,y_hat,'bo'); axis([-5,5,-5,5]); hold on
plot([-5, 5], [-5,5])
saveas(gcf,'fig2b.jpg')

%% shuffle and split crossvalidation
nfolds = 10000;
training_mads = zeros(nfolds,1);
test_mads = zeros(nfolds,1);
training_rmses = zeros(nfolds,1);
test_rmses = zeros(nfolds,1);
adj_rs = zeros(nfolds,1);
training_resids = zeros(nfolds,1258);
training_inds = zeros(nfolds,1258);

test_resids = zeros(nfolds,139);
test_inds = zeros(nfolds,139);

for fold = 1:nfolds
	fold
	reg_table.inds = (1:length(reg_table.distance))';
	reg_table.rand = randperm(length(reg_table.distance))';
	reg_table = sortrows(reg_table, width(reg_table));
	reg_table = removevars(reg_table,'rand');
	test_inds(fold,:)= reg_table.inds(1:139);
	training_inds(fold,:)= reg_table.inds(140:end);
	reg_table = removevars(reg_table,'inds');
	
	test_data = reg_table(1:139,:);
	training_data = reg_table(140:end,:);
	training_mdl = fitlm(training_data,'interactions','ResponseVar', 'distance', 'RobustOpts','on');
	training_resids(fold,:) = training_mdl.Residuals.Raw;
	test_resids(fold,:) = test_data.distance-predict(training_mdl, test_data);
	adj_rs(fold) = training_mdl.Rsquared.Adjusted;
end
training_mads = nanmedian(abs(training_resids),2);
test_mads = nanmedian(abs(test_resids),2);
mad_ratios = training_mads./test_mads;
save('cross_val_results.mat','adj_rs', 'training_resids','test_resids', 'training_mads',...
	'test_mads', 'mad_ratios');

%% cross-validate figures
load cross_val_results
figure; % figure 3
hist(adj_rs,100);
mean(adj_rs)
median(adj_rs)
mean(adj_rs(adj_rs>.8))
saveas(gcf,'fig3.jpg')

figure; % figure 4
subplot(3,1,1); hist(training_mads,40); xlim([0 1.5])
subplot(3,1,2); hist(test_mads,40); xlim([0 1.5])
subplot(3,1,3); hist(mad_ratios,40); xlim([0 1.5])
saveas(gcf,'fig4.jpg')

sum(mdl_robust.Coefficients.pValue<0.05)
