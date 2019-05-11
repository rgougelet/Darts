clear; close all; clc;
script_dir = 'C:\Users\Rob\Desktop\darts\';
cd(script_dir);
addpath('.\eeglab13_6_5b')
data_dir = '.\data\';
addpath(data_dir);

% user input
srate = 64;
chan_labs = {'Fz','Cz','Pz','Oz'};
freq_labs = {'theta', 'alpha'};
interval_labs = {'encoding','maintenance', 'recall'};
subjs_to_include = {'571', '579', '580', ...
	'607', '608', '616', '619', '621', '627', '631'};

% load XLSX SNAP data
[num,txt,raw] = xlsread([data_dir,'behavioral_data.xlsx']);
headers = txt(1,:);
for k=1:numel(headers)
	xlsx.(headers{k})=num(:,k) ;
end

% initialize eeg columns to add to the behav columns
for chan_lab = chan_labs
	for freq_lab = freq_labs
		for interval_lab = interval_labs
			xlsx.([interval_lab{:},'_',chan_lab{:},'_',freq_lab{:}]) = nan(length(xlsx.pres),1);
		end
	end
end

for subj_i = 1:length(subjs_to_include)
	subj_id = subjs_to_include{subj_i};
	subj_set = [subj_id,'_eeg_64_laplace.set'];
	
	% bandpass filter at desired freqs
	EEG = pop_loadset('filename',subj_set,'filepath',data_dir);
	EEG_theta = pop_eegfiltnew(EEG, 3, 7);
	EEG_alpha = pop_eegfiltnew(EEG, 8, 12);
	EEGs = {EEG_theta,EEG_alpha};
	load([data_dir, subj_id,'_eeg_64_latencys'])
	
	% correct data collection issues
	% the problem is there are some trials in the
	% xls file that are not in eeg
	n_snap_trials = sum(xlsx.subject == str2double(subj_id));
	n_eeg_trials = length(end_event_latencys);
	eeg_trial_strs = str2num(end_event_strings(:,1:4)); % ignore warning, use str2num
	subj_inds = xlsx.subject == str2double(subj_id);
	snap_trial_strs = str2num([num2str(xlsx.delay(subj_inds)),...
		num2str(xlsx.position(subj_inds),'%02d'),num2str(xlsx.pres(subj_inds))]);
	eeg_to_snap_inds = 1:length(eeg_trial_strs);
	if strcmp(subj_id, '580')
		eeg_to_snap_inds = 10+(1:n_eeg_trials);
	end
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

	% retrieve this subjs trial-level amplitudes at desired chans and freqs
	for chan_lab_i = 1:length(chan_labs)
		chan_lab = chan_labs{chan_lab_i};
		for freq_lab_i = 1:length(EEGs)
			freq_lab = freq_labs{freq_lab_i};
			EEG = EEGs{freq_lab_i};
% 			encode_trial_amps = [];
			maint_trial_amps = [];
			recall_trial_amps = [];
			for eeg_trial_i = 1:n_eeg_trials
				start_latency_i = start_event_latencys(eeg_trial_i);
				cue_latency_i = cue_event_latencys(eeg_trial_i);
				end_latency_i = end_event_latencys(eeg_trial_i);
% 				encode_trial_amps(end+1) = mean(abs(EEG.data(chan_lab_i,start_latency_i:(start_latency_i+srate))),2)';
				maint_baseline = mean((EEG.data(chan_lab_i,round(start_latency_i-0.2*srate):start_latency_i)).^2,2)';
				maint_trial_amps(end+1) = sum((EEG.data(chan_lab_i,start_latency_i:cue_latency_i)).^2,2)'-maint_baseline;
				recall_baseline = mean((EEG.data(chan_lab_i,round(cue_latency_i-0.2*srate):cue_latency_i)).^2,2)';
				recall_trial_amps(end+1) = sum((EEG.data(chan_lab_i,cue_latency_i:end_latency_i)).^2,2)'-recall_baseline;
			end
			
% 			encode_chan_freq_amps = xlsx.(['encoding_',chan_lab,'_',freq_lab]);
			maint_amps = xlsx.(['maintenance_',chan_lab,'_',freq_lab]);
			recall_amps = xlsx.(['recall_',chan_lab,'_',freq_lab]);

			for trial_amp_i = 1:n_eeg_trials
% 				encode_chan_freq_amps(eeg_to_snap_inds(trial_amp_i)) = encode_trial_amps(trial_amp_i);
				maint_amps(eeg_to_snap_inds(trial_amp_i)) = maint_trial_amps(trial_amp_i);
				recall_amps(eeg_to_snap_inds(trial_amp_i)) = recall_trial_amps(trial_amp_i);
			end
% 			xlsx.(['encoding_',chan_lab,'_',freq_lab]) = encode_chan_freq_amps;
			xlsx.(['maintenance_',chan_lab,'_',freq_lab]) = maint_amps;
			xlsx.(['recall_',chan_lab,'_',freq_lab]) = recall_amps;
	
		end
	end	
end

%%
% parsave('behav+eeg.mat', xlsx, 'xlsx')

%% behav+eeg regression, no interaction 
non_outlier_inds = (xlsx.distance < 6 & xlsx.throwtime < 5);
throwtime = xlsx.throwtime(non_outlier_inds);
delay = xlsx.delaystilltime(non_outlier_inds);
pres = (xlsx.pres(non_outlier_inds));
dists = xlsx.distance(non_outlier_inds).^(1/2);
subjs = categorical(xlsx.subject(non_outlier_inds));


chan_labs = {'Fz','Cz','Pz','Oz'};
freq_labs = {'theta', 'alpha'};
interval_labs = {'maintenance', 'recall'};
varnames = {'dists'; 'throwtime'; 'delay';'pres'};
for chan_lab = chan_labs
	for freq_lab = freq_labs
		for interval_lab = interval_labs
			varname = [interval_lab{:},'_',chan_lab{:},'_',freq_lab{:}];
			varnames{end+1} = varname;
			eval([varname,' = xlsx.',varname,'(non_outlier_inds);']);
		end
	end
end

eval( ['array = [', strjoin(varnames,','),'];'] );
[rho,pval] = corr(array, 'rows','complete');
heatmap(rho)
heatmap(pval)
%%
eval( ['reg_table = table(', strjoin(varnames,','),');'] );

fit_eeg = fitlm(reg_table,'interactions','ResponseVar', 'dists', 'RobustOpts', 'on')
results = sortrows(anova(fit_eeg),5,'descend')
% results = table2array(anova(fit_eeg));
% results(results(:,5) < 0.05)
fit_eeg = stepwiselm(reg_table,'interactions','ResponseVar', 'dists')
% 
% thresholds = 0.05./(height(results):-1:1)'
% sum(results.pValue<thresholds)
% fdr_res = fdr_bh(results.pValue, 0.9999, 'dep', 'yes');
% for k = 1:length(results.pValue)
% 	if results.pValue(k) <= (k/length(results.pValue))*0.05
% 		break
% 	end
% end
% m = length(results.pValue);
% results.pValue<=0.05
% plot(results.pValue)
%%
for varname = varnames(4:end)
	varnames{end+1} = [varname{1};'*throwtime']
	varnames{end+1} = [varname{1},'*delay']
end

fit_eeg = fitlm(reg_table, ['dists ~ ',strjoin(varnames(2:end),' + ')], 'RobustOpts', 'on')
fit_eeg = stepwiselm(reg_table, ['dists ~ ',strjoin(varnames(2:end),' + ')])


% % [H,P,CI,STATS] = ttest2(Cz_thetas_abs,Cz_thetas_pres)
% % [H,P,CI,STATS] = ttest2(Cz_thetas_abs,Cz_thetas_pres, 'vartype', 'unequal')

%% behavioral regression
varnames = {'dists', 'throwtime', 'delay', 'pres'};
eval( ['reg_table = table(', strjoin(varnames,','),');'] );
fit_behav = fitlm(reg_table,'interactions', 'ResponseVar', 'dists', 'RobustOpts', 'on')
swfit_behav = stepwiselm(reg_table,'interactions', 'ResponseVar', 'dists', 'Criterion', 'AdjRsquared')
% anova(swfit_behav)