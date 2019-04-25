clear; close all; clc;
script_dir = 'C:\Users\Rob\Desktop\darts\';
cd(script_dir);
addpath('.\eeglab13_6_5b')
data_dir = '.\data\';
addpath(data_dir);
eeglab;
close all;

% user input
chan_labs = {'Fz','Cz','Pz','Oz'};
freq_labs = {'theta', 'alpha'};
% subjs_to_include = {'571', '579', '580', ...
% 	'607', '608', '616', '619', '621', '627', '631'};
subjs_to_include = {'616'};


% load XLSX SNAP data
[num,txt,raw] = xlsread([data_dir,'behavioral_data.xlsx']);
headers = txt(1,:);
for k=1:numel(headers)
	xlsx.(headers{k})=num(:,k) ;
end

% initialize eeg columns to add the behav columns
for chan_lab = chan_labs
	for freq_lab = freq_labs
		xlsx.([chan_lab{:},'_',freq_lab{:}]) = nan(length(xlsx.pres),1);
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
	eeg_trial_strs = str2num(end_event_strings(:,1:4)); % ignore warning, use str2num
	n_eeg_trials = length(end_event_time_inds);
	subj_inds = xlsx.subject == str2double(subj_id);
	snap_trial_strs = str2num([num2str(xlsx.delay(subj_inds)),num2str(xlsx.position(subj_inds),'%02d'),num2str(xlsx.pres(subj_inds))]);
	n_snap_trials = length(snap_trial_strs);
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

	% retrieve this subjs trial-level amplitudes at desired chans and freqs
	for chan_lab_i = 1:length(chan_labs)
		chan_lab = chan_labs{chan_lab_i};
		for freq_lab_i = 1:length(EEGs)
			freq_lab = freq_labs{freq_lab_i};
			EEG = EEGs{freq_lab_i};
			trial_amps = [];
			for eeg_trial_i = 1:length(start_event_time_inds)
				start_latency_i = start_event_time_inds(eeg_trial_i);
				end_latency_i = end_event_time_inds(eeg_trial_i);
				trial_amps(end+1) = mean(abs(EEG.data(chan_lab_i,start_latency_i:end_latency_i)),2)';
			end
			all_amps = xlsx.([chan_lab,'_',freq_lab]);
			all_amps(xl_start_row:(xl_start_row+length(start_event_time_inds)-1)) = trial_amps;
			xlsx.([chan_lab,'_',freq_lab]) = all_amps;
		end
	end	
end

% trial-level behav+eeg regression
throwtimes = xlsx.throwtime(xlsx.distance < 6 & xlsx.throwtime < 5);
delay = xlsx.delaystilltime(xlsx.distance < 6 & xlsx.throwtime < 5);
pres = categorical(xlsx.pres(xlsx.distance < 6 & xlsx.throwtime < 5));
dists = xlsx.distance(xlsx.distance < 6 & xlsx.throwtime < 5);
subjs = categorical(xlsx.subject(xlsx.distance < 6 & xlsx.throwtime < 5));

% create chan+freq  variables
varnames = {'dists', 'throwtimes', 'delay', 'pres', 'subjs'};
for chan_lab = chan_labs
	for freq_lab = freq_labs
		varname = [chan_lab{:},'_',freq_lab{:}];
		varnames{end+1} = varname;
		eval([varname,' = xlsx.',varname,'(xlsx.distance < 6 & xlsx.throwtime < 5);']);
	end
end

eval( ['reg_table = table(', strjoin(varnames,','),');'] );

fit = fitlm(reg_table,'ResponseVar', 'dists', 'RobustOpts', 'on');


% [H,P,CI,STATS] = ttest2(Cz_thetas_abs,Cz_thetas_pres)
% [H,P,CI,STATS] = ttest2(Cz_thetas_abs,Cz_thetas_pres, 'vartype', 'unequal')